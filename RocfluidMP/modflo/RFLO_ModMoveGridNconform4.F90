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
! Purpose: Suite for frame grid-motion routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModMoveGridNconform4.F90,v 1.13 2009/08/27 14:04:51 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModMoveGridFrame

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
  PUBLIC :: RFLO_MoveGridFrame, &
            RFLO_MgFrameCornPoints, &
            RFLO_MgFrameBroadcast, &
            RFLO_MgFrameSrchNeighbors, &

            RFLO_MgFrameCorrectNeighbors, &
            RFLO_MgFrameMoveCorners, &
            RFLO_MgFrameOrthoShift, &
            RFLO_MgFrameEdges, &
            RFLO_MgFrameRestoreExternal, &
            RFLO_MgFrameBndDeformation

! private : RFLO_MgFrameSurface
!           RFLO_MgFrameInterfaces
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModMoveGridNconform4.F90,v $ $Revision: 1.13 $'        
             
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

SUBROUTINE RFLO_MoveGridFrame( regions )

  USE ModInterfaces, ONLY : RFLO_MoveGridSurfaces, RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, RFLO_LaplaceGridSmoo,&
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
        RFLO_CheckMetrics, RFLO_CalcGridSpeeds, RFLO_BoundaryDeformation

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

  INTEGER :: bcType, iRemesh, jRemesh, nRemesh, iType

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

  CALL RegisterFunction( global,'RFLO_MoveGridFrame',&
  'RFLO_ModMoveGridNconform4.F90' )

  iType=1

#ifdef GENX
! update geometry buffers -----------------------------------------------------

  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function( global%genxHandleGm,1,dAlpha )
#endif

! receive and distribute deformations for each region -------------------------

  CALL RFLO_MgFrameSurfaces( regions,someMoved,iType )

! fix interfaces between regions ----------------------------------------------

  IF (someMoved) THEN
    CALL RFLO_MgFrameInterfaces( regions,iType )
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

! smooth grid by solving Laplace equation -------------------------------------

  IF (global%moveGridNiter < 1) THEN 
    IF (global%verbLevel >= VERBOSE_HIGH) THEN
      IF (global%myProcid == MASTERPROC) THEN
        WRITE(STDOUT,3000) SOLVER_NAME,global%skewness,global%minVol
        WRITE(STDOUT,1000) SOLVER_NAME, &
                  global%moveGridNiter,global%moveGridNsmatch, &
                  global%moveGridAmplifX,global%moveGridAmplifY, &
                  global%moveGridAmplifZ,global%moveGridPower, &
                  global%moveGridOrthDir,global%moveGridOrthWghtX, &
                  global%moveGridOrthWghtY,global%moveGridOrthWghtZ
      ENDIF  ! masterproc
    ENDIF    ! verbLevel
    GOTO 888
  ENDIF      ! niter<1

  IF (someMoved) THEN
    DO iter=1,global%moveGridNiter
      CALL RFLO_LaplaceGridSmoo( regions,resid )
    ENDDO

    IF (global%verbLevel >= VERBOSE_HIGH) THEN
#ifdef MPI
      CALL MPI_Reduce( resid,globalResid,1,MPI_RFREAL,MPI_SUM, &
                       MASTERPROC,global%mpiComm,global%mpierr )
      IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#else
      globalResid = resid
#endif
      IF (global%myProcid == MASTERPROC) THEN
        WRITE(STDOUT,3000) SOLVER_NAME,global%skewness,global%minVol
        WRITE(STDOUT,2000) SOLVER_NAME, &
                  global%moveGridNiter,global%moveGridNsmatch, &
                  global%moveGridAmplifX,global%moveGridAmplifY, &
                  global%moveGridAmplifZ,global%moveGridPower, &
                  global%moveGridOrthDir,global%moveGridOrthWghtX, &
                  global%moveGridOrthWghtY,global%moveGridOrthWghtZ, &
                  SQRT(globalResid)
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
        IF ((bcType==BC_SYMMETRY)) THEN
           grid%boundMoved(patch%lbound) = .FALSE.
        ENDIF  ! bcType
     ENDDO    ! iPatch
!      CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
!                                     grid%edgeMoved,grid%arcLen12, &
!                                     grid%arcLen34,grid%arcLen56, &
!                                     gridOld%xyzOld,grid%xyz )

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

888 CONTINUE

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

! global grid quality measure -------------------------------------------------

  CALL RFLO_GridQualityGlobal( regions )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(A,1X,'Global-Weighted-Laplacian grid motion:', &
            2I5,4(1PE9.2),I4,3(1PE9.2))
2000 FORMAT(A,1X,'Global-Weighted-Laplacian grid motion:', &
            2I5,4(1PE9.2),I4,3(1PE9.2),1PE13.4)
3000 FORMAT(A,1X,'global skewness, minvol:',2(1PE14.5))

END SUBROUTINE RFLO_MoveGridFrame

!******************************************************************************
!
! Purpose: search for corner points including those of internal patches
!
! Description: none.
!
! Input: regions = data of current region.
!
! Output: grid%nCorns  = number of corner points in each region
!         grid%ijkCorn = ijkValue of each corner
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameCornPoints( regions )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, &
                            RFLO_GetDimensPhysNodes, RFLO_GetPatchIndicesNodes

  IMPLICIT NONE
#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: l, iPatch, iReg, ipCorn, intCorn, nReg

! ... local variables
  INTEGER, PARAMETER :: ncMax=100
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: iptc, jptc, kptc, iblk, jblk, kblk, ijkCurr
  INTEGER :: iNOff, ijNOff, lbound, regNc, errFl
  INTEGER, ALLOCATABLE :: ivar(:), ijkCorn(:,:)
  LOGICAL :: wasfound

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameCornPoints',&
  'RFLO_ModMoveGridNconform4.F90' )

! search for block and patch corners in each region ---------------------------

  ALLOCATE( ivar(global%nRegions),stat=errFl )
  global%error = errFl
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( ijkCorn(ncMax,global%nRegions),stat=errFl )
  global%error = errFl
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  iLev = 1

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

! --- search for internal patch corners

      grid%nCorns(iReg) = 8
      ijkCorn(1,iReg) = IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff)
      ijkCorn(2,iReg) = IndIJK(ipnbeg,jpnbeg,kpnend,iNOff,ijNOff)
      ijkCorn(3,iReg) = IndIJK(ipnbeg,jpnend,kpnend,iNOff,ijNOff)
      ijkCorn(4,iReg) = IndIJK(ipnbeg,jpnend,kpnbeg,iNOff,ijNOff)
      ijkCorn(5,iReg) = IndIJK(ipnend,jpnbeg,kpnbeg,iNOff,ijNOff)
      ijkCorn(6,iReg) = IndIJK(ipnend,jpnbeg,kpnend,iNOff,ijNOff)
      ijkCorn(7,iReg) = IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff)
      ijkCorn(8,iReg) = IndIJK(ipnend,jpnend,kpnbeg,iNOff,ijNOff)

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(iLev)%patches(iPatch)
        lbound =  patch%lbound

        CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                        ibeg,iend,jbeg,jend,kbeg,kend )

        DO ipCorn = 1,4    ! patch corners
          IF (lbound==1 .OR. lbound==2) THEN
            iptc = ibeg
            IF (lbound==1) iblk = ipnbeg
            IF (lbound==2) iblk = ipnend
            IF (ipCorn==1) THEN
              jptc = jbeg
              jblk = jpnbeg
              kptc = kbeg
              kblk = kpnbeg
            ELSEIF (ipCorn==2) THEN
              jptc = jbeg
              jblk = jpnbeg
              kptc = kend
              kblk = kpnend
            ELSEIF (ipCorn==3) THEN
              jptc = jend
              jblk = jpnend
              kptc = kend
              kblk = kpnend
            ELSEIF (ipCorn==4) THEN
              jptc = jend
              jblk = jpnend
              kptc = kbeg
              kblk = kpnbeg
            ENDIF
          ELSEIF (lbound==3 .OR. lbound==4) THEN
            jptc = jbeg
            IF (lbound==3) jblk = jpnbeg
            IF (lbound==4) jblk = jpnend
            IF (ipCorn==1) THEN
              kptc = kbeg
              kblk = kpnbeg
              iptc = ibeg
              iblk = ipnbeg
            ELSEIF (ipCorn==2) THEN
              kptc = kbeg
              kblk = kpnbeg
              iptc = iend
              iblk = ipnend
            ELSEIF (ipCorn==3) THEN
              kptc = kend
              kblk = kpnend
              iptc = iend
              iblk = ipnend
            ELSEIF (ipCorn==4) THEN
              kptc = kend
              kblk = kpnend
              iptc = ibeg
              iblk = ipnbeg
            ENDIF  ! ipCorn
          ELSEIF (lbound==5 .OR. lbound==6) THEN
            kptc = kbeg
            IF (lbound==5) kblk = kpnbeg
            IF (lbound==6) kblk = kpnend
            IF (ipCorn==1) THEN
              iptc = ibeg
              iblk = ipnbeg
              jptc = jbeg
              jblk = jpnbeg
            ELSEIF (ipCorn==2) THEN
              iptc = ibeg
              iblk = ipnbeg
              jptc = jend
              jblk = jpnend
            ELSEIF (ipCorn==3) THEN
              iptc = iend
              iblk = ipnend
              jptc = jend
              jblk = jpnend
            ELSEIF (ipCorn==4) THEN
              iptc = iend
              iblk = ipnend
              jptc = jbeg
              jblk = jpnbeg
            ENDIF  ! ipCorn
          ENDIF    ! lbound

          patch%corns(ipCorn) = IndIJK(iptc,jptc,kptc,iNOff,ijNOff) 

          IF (iptc/=iblk .OR. jptc/=jblk .OR. kptc/=kblk) THEN
            wasfound = .false.
            ijkCurr  = IndIJK(iptc,jptc,kptc,iNOff,ijNOff)
            DO intCorn=1,grid%nCorns(iReg)
              IF (ijkCorn(intCorn,iReg)==ijkCurr) THEN
                wasfound = .true.
              ENDIF
            ENDDO
            IF (.NOT. wasfound) THEN
              grid%nCorns(iReg) = grid%nCorns(iReg) +1
              ijkCorn(grid%nCorns(iReg),iReg) = ijkCurr
            ENDIF
          ENDIF
          IF (grid%nCorns(iReg) >= ncMax) THEN
            CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
            'too low ncMax in RFLO_ModMoveGridFrame/RFLO_MgFrameCornPoints')
          ENDIF
        ENDDO   ! ipCorn
      ENDDO   ! iPatch

      ivar(iReg) = grid%nCorns(iReg)
    ENDIF  ! myProcid
  ENDDO    ! iReg

#ifdef MPI
  DO iReg = 1,global%nRegions
    CALL MPI_BCAST( ivar(iReg),1,MPI_INTEGER, &
                    regions(iReg)%procId,global%mpiComm,global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

    CALL MPI_BCAST( ijkCorn(1:ncMax,iReg),ncMax,MPI_INTEGER, &
                    regions(iReg)%procId,global%mpiComm,global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  ENDDO
#endif

  regNc = 0
  DO iReg = 1,global%nRegions
    regNc = MAX( regNc,ivar(iReg) )
  ENDDO
  global%moveGridRegNc = regNc

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid
      DO nReg = 1,global%nRegions
        grid%nCorns(nReg) = ivar(nReg)
      ENDDO

      ALLOCATE( grid%ijkCorn(      regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%regCorn(    3,regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%regCornOld( 3,regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%regCornOrig(3,regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%nghbor(     3,global%moveGridNbour,regNc),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      DO l = 1,grid%nCorns(iReg)
        grid%ijkCorn(l,iReg) = ijkCorn(l,iReg)
      ENDDO

! --- shared corners

      ALLOCATE( grid%nShared( regNc),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%regCornBuff(3,regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ALLOCATE( grid%regCornOrth(3,regNc,global%nRegions),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ENDIF  ! myProcid
  ENDDO    ! iReg

! deallocate temporary arrays

  DEALLOCATE( ivar,stat=errFl )
  global%error = errFl
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameCornPoints

!******************************************************************************
!
! Purpose: broadcast movements at 8 corner points of current region to all 
!          regions
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Notes: upon first call by RFLO_InitGridProcedure, regions%levels%grid%xyz 
!        contains grid coordinates, but on second call by RFLO_MgFrameSurface
!        regions%levels%grid%xyz contains grid movements.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameBroadcast( regions,iselect,iter )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: iselect, iter

! ... loop variables
  INTEGER :: i, l, iReg

! ... local variables
  INTEGER :: iLev, nCorns, errFl
  INTEGER, ALLOCATABLE :: corner(:)

  REAL(RFREAL), ALLOCATABLE :: rvar(:,:,:)
  REAL(RFREAL), POINTER :: dxyz(:,:), xyzOld(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER   :: grid, gridOld

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameBroadcast',&
  'RFLO_ModMoveGridNconform4.F90' )

! store block corners and broadcast to all regions ----------------------------

  ALLOCATE( rvar(XCOORD:ZCOORD,global%moveGridRegNc,global%nRegions), &
            stat=errFl )
  global%error = errFl
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  rvar = 0._RFREAL

  ilev = 1

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid    => regions(iReg)%levels(iLev)%grid
      gridOld => regions(iReg)%levels(iLev)%gridOld

      dxyz   => grid%xyz
      xyzOld => gridOld%xyz
      nCorns =  grid%nCorns(iReg)

      ALLOCATE( corner(nCorns), stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      DO l = 1,nCorns
        corner(l) = grid%ijkCorn(l,iReg)
      ENDDO

      IF (iter==1) THEN
!        CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
!                                      jpnbeg,jpnend,kpnbeg,kpnend )
!        CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

        IF (iselect==0) THEN
          DO i=1,nCorns
            grid%regCornOrig(XCOORD,i,iReg) = dxyz(XCOORD,corner(i))
            grid%regCornOrig(YCOORD,i,iReg) = dxyz(YCOORD,corner(i))
            grid%regCornOrig(ZCOORD,i,iReg) = dxyz(ZCOORD,corner(i))
            rvar(:,i,iReg) = grid%regCornOrig(:,i,iReg)
          ENDDO
        ELSEIF (iselect==1) THEN
          DO i=1,nCorns
            grid%regCornOld(XCOORD,i,iReg) = dxyz(XCOORD,corner(i))
            grid%regCornOld(YCOORD,i,iReg) = dxyz(YCOORD,corner(i))
            grid%regCornOld(ZCOORD,i,iReg) = dxyz(ZCOORD,corner(i))
            rvar(:,i,iReg) = grid%regCornOld(:,i,iReg)
          ENDDO
        ENDIF

      ELSE
        IF (iselect==1) &
        rvar(:,:,iReg) = grid%regCornOld(:,:,iReg)

      ENDIF  ! iter

      IF (iselect==2) THEN
        DO i=1,nCorns
!          grid%regCornBuff(:,i,iReg) = xyzOld(:,corner(i))+dxyz(:,corner(i))
          grid%regCornBuff(:,i,iReg) = xyzOld(:,corner(i))
          rvar(:,i,iReg) = grid%regCornBuff(:,i,iReg)
        ENDDO
      ENDIF
      IF (iselect==3) THEN
        DO i=1,nCorns
          rvar(:,i,iReg) = grid%regCornOrth(:,i,iReg)
        ENDDO
      ENDIF

      DEALLOCATE( corner, stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

    ENDIF    ! myProcid
  ENDDO      ! iReg

#ifdef MPI
  DO iReg = 1,global%nRegions
    nCorns =  global%moveGridRegNc

    CALL MPI_BCAST( rvar(XCOORD:ZCOORD,1:nCorns,iReg),3*nCorns, &
               MPI_RFREAL,regions(iReg)%procId,global%mpiComm,global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
  ENDDO
  CALL MPI_Barrier( global%mpiComm,global%mpierr )

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid
      IF (iter == 1) THEN
        IF (iselect==0) THEN
          DO l=1,global%nRegions
            grid%regCornOrig(:,:,l) = rvar(:,:,l)
          ENDDO
        ELSEIF (iselect==1) THEN
          DO l=1,global%nRegions
            grid%regCornOld(:,:,l) = rvar(:,:,l)
          ENDDO
        ENDIF
      ELSE
        IF (iselect==1) THEN
          DO l=1,global%nRegions
            grid%regCornOld(:,:,l) = rvar(:,:,l)
          ENDDO
        ENDIF
      ENDIF  ! iter
      IF (iselect==2) THEN
        DO l=1,global%nRegions
          grid%regCornBuff(:,:,l) = rvar(:,:,l)
        ENDDO
      ENDIF
      IF (iselect==3) THEN
        DO l=1,global%nRegions
          grid%regCornOrth(:,:,l) = rvar(:,:,l)
        ENDDO
      ENDIF  ! iselect
    ENDIF    ! myProcid
  ENDDO      ! iReg

!  DO iReg = 1,global%nRegions
!    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
!        regions(iReg)%active==ACTIVE) THEN               ! on my processor
!      grid => regions(1)%levels(iLev)%grid
!      DO l = 1,global%nRegions
!        DO i=1,grid%nCorns(iReg)
!          write(*,*)iReg,l,i,grid%regCornOrig(:,i,l)
!        ENDDO
!      ENDDO
!    ENDIF
!  ENDDO
#endif

! deallocate temporary arrays

  DEALLOCATE( rvar,stat=errFl )
  global%error = errFl
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameBroadcast

!******************************************************************************
!
! Purpose: search for six closest neighbors
!
! Description: none.
!
! Input: regions = data of current region.
!
! Output: grid%nghbor = neighbouring points identified
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameSrchNeighbors( regions )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetNodeOffset, &
                            RFLO_GetDimensPhysNodes, RFLO_GetPatchIndicesNodes

  IMPLICIT NONE
#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: i, j, k, iPatch, ic, iReg, nc, nReg

! ... local variables
  INTEGER, PARAMETER :: nSharedMax=16
  INTEGER :: iLev, bcType, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: regNc, nCorns, nbour, nshared, nsmax
  INTEGER :: ijkNode(4), iNOff, ijNOff, lbound, errFl
  INTEGER, ALLOCATABLE :: ncMin(:), nRegMin(:), cshared(:,:,:)
  REAL(RFREAL) :: edgeLen, ds, tol
  REAL(RFREAL), POINTER :: xyz(:,:)
  REAL(RFREAL), ALLOCATABLE :: dist(:,:), distMin(:)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameSrchNeighbors',&
  'RFLO_ModMoveGridNconform4.F90' )

! search for six closest neighbours -------------------------------------------

  global%MoveGridNsharedMax = nSharedMax  ! max number of shared neigbors
  nbour = global%moveGridNbour            ! number of closest neighbors
  regNc = global%moveGridRegNc            ! max number of block-corners 
  iLev  = 1

  ALLOCATE( ncMin(nbour), stat=errFl );                IF (errFl>0) GOTO 88
  ALLOCATE( nRegMin(nbour), stat=errFl );              IF (errFl>0) GOTO 88
  ALLOCATE( distMin(nbour), stat=errFl );              IF (errFl>0) GOTO 88
  ALLOCATE( dist(regNc,global%nRegions), stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( cshared(2,nSharedMax,regNc), stat=errFl ); IF (errFl>0) GOTO 88

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      nCorns = grid%nCorns(iReg)

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

      xyz => regions(iReg)%levels(iLev)%grid%xyz

! --- calculate the shortest cell edge

      edgeLen = 1.E+30_RFREAL

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ijkNode(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
            ijkNode(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
            ijkNode(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
            ijkNode(4) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
            ds      = SQRT((xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
            ds      = SQRT((xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
            ds      = SQRT((xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
          ENDDO
        ENDDO
      ENDDO
      tol = 1.E-5_RFREAL*edgeLen

      nsmax = 0

      DO ic = 1,nCorns
        distMin(1:nbour) = 1.e+30_RFREAL
        ncMin(1:nbour)   = 1
        nRegMin(1:nbour) = 1
        nshared          = 0
        DO nReg = 1,global%nRegions
          DO nc = 1,grid%nCorns(nReg)
            dist(nc,nReg) = SQRT((grid%regCornOrig(XCOORD,nc,nReg)- &
                                  grid%regCornOrig(XCOORD,ic,iReg))**2 + &
                                 (grid%regCornOrig(YCOORD,nc,nReg)- &
                                  grid%regCornOrig(YCOORD,ic,iReg))**2 + &
                                 (grid%regCornOrig(ZCOORD,nc,nReg)- &
                                  grid%regCornOrig(ZCOORD,ic,iReg))**2)

! -------- inhibitor check
!           IF (dist(nc,nReg)>edgeLen .AND. iReg==12 .AND. ic==5 .AND. &
!               (nReg==12 .OR. nReg==14 .OR. nReg==18 .OR. nReg==25 .OR. &
!                nReg==26 .OR. nReg==52 .OR. nReg==69)) &
!             write(*,*)'i',iReg,ic,nReg,nc,dist(nc,nReg)

! -------- titan4-240blocks check
!           IF (dist(nc,nReg)>edgeLen .AND. iReg==120 .AND. ic==6 .AND. &
!               (nReg==96 .OR. nReg==37 .OR. nReg==105 .OR. nReg==121 .OR. &
!                nReg==117 .OR. nReg==123 .OR. nReg==71 .OR. nReg==72 .OR. &
!                nReg==120)) write(*,*)'i',iReg,ic,nReg,nc,dist(nc,nReg)

            IF (dist(nc,nReg)<distMin(1) .AND. dist(nc,nReg)>edgeLen) THEN
              DO k = nbour,2,-1
                distMin(k) = distMin(k-1)
                ncMin(k)   = ncMin(k-1)
                nRegMin(k) = nRegMin(k-1)
              ENDDO
              distMin(1) = dist(nc,nReg)
              ncMin(1)   = nc
              nRegMin(1) = nReg
            ENDIF

            DO k = 2,nbour
              IF (dist(nc,nReg) > (distMin(k-1) + tol) .AND. &
                  dist(nc,nReg) < (distMin(k) - tol)) THEN
! ------------- titan4-240blocks check
!               IF (iReg==120 .AND. ic==6) write(*,*)'j',iReg,ic,nReg, &
!                 nc,k,dist(nc,nReg),distMin(k-1),distMin(k)
                DO j = nbour,k+1,-1
                  distMin(j) = distMin(j-1)
                  ncMin(j)   = ncMin(j-1)
                  nRegMin(j) = nRegMin(j-1)
                ENDDO
                distMin(k) = dist(nc,nReg)
                ncMin(k)   = nc
                nRegMin(k) = nReg
              ENDIF
            ENDDO

! --------- search for neigboring block/corners shared with this corner
            IF (dist(nc,nReg)<edgeLen) THEN
              nshared = nshared+1
              IF (nshared > nSharedMax) THEN
                CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
                'too low nSharedMax in RFLO_MgFrameSrchNeighbors')
              ENDIF
              cshared(1,nshared,ic) = nc
              cshared(2,nshared,ic) = nReg
            ENDIF

          ENDDO  ! nc
        ENDDO    ! nReg

        grid%nShared(ic) = nshared
        nsmax = MAX( nsmax,nshared )

! ----- store neighbor info

        DO k = 1,nbour
          nc   = ncMin(k)
          nReg = nRegMin(k)
          grid%nghbor(1,k,ic) = nc
          grid%nghbor(2,k,ic) = nReg
        ENDDO    ! k
      ENDDO      ! ic

      grid%nghbor(3,:,:) = 1

! --- store shared corners info of current region

      ALLOCATE( grid%cshared(2,nsmax,regNc),stat=errFl )
      global%error = errFl
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      grid%cshared(1:2,1:nsmax,1:nCorns) = cshared(1:2,1:nsmax,1:nCorns)

    ENDIF  ! myProcid
  ENDDO    ! iReg

! deallocate temporary arrays -------------------------------------------------

  DEALLOCATE( ncMin,   stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( nRegMin, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( distMin, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( dist,    stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( cshared, stat=errFl ); IF (errFl>0) GOTO 99

! assign internal/external flag to block corners

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

! --- corner 1
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(1)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(1)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(4)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(9)%interact)) &
        grid%nghbor(3,1:nbour,1) = 0

! --- corner 2
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(2)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(1)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(2)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(10)%interact)) &
        grid%nghbor(3,1:nbour,2) = 0

! --- corner 3
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(3)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(2)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(3)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(11)%interact)) &
        grid%nghbor(3,1:nbour,3) = 0

! --- corner 4
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(4)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(3)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(4)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(12)%interact)) &
        grid%nghbor(3,1:nbour,4) = 0

! --- corner 5
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(5)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(5)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(8)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(9)%interact)) &
        grid%nghbor(3,1:nbour,5) = 0

! --- corner 6
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(6)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(5)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(6)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(10)%interact)) &
        grid%nghbor(3,1:nbour,6) = 0

! --- corner 7
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(7)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(6)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(7)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(11)%interact)) &
        grid%nghbor(3,1:nbour,7) = 0

! --- corner 8
      IF ((.NOT. regions(iReg)%levels(iLev)%cornerCells(8)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(7)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(8)%interact).OR. &
          (.NOT. regions(iReg)%levels(iLev)%edgeCells(12)%interact)) &
        grid%nghbor(3,1:nbour,8) = 0

    ENDIF  ! myProcid
  ENDDO    ! iReg

! second step search for external corners

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(iLev)%patches(iPatch)
        lbound =  patch%lbound
        bcType =  patch%bcType

        CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                        ibeg,iend,jbeg,jend,kbeg,kend )

        IF ((bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
            (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR. &
            (bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR. &
            (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
            (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
            (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          IF (lbound==1 .OR. lbound==2) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)
            ijkNode(3) = IndIJK(ibeg,jend,kend,iNOff,ijNOff)
            ijkNode(4) = IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)
          ELSEIF (lbound==3 .OR. lbound==4) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)
            ijkNode(3) = IndIJK(iend,jbeg,kend,iNOff,ijNOff)
            ijkNode(4) = IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)
          ELSEIF (lbound==5 .OR. lbound==6) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)
            ijkNode(3) = IndIJK(iend,jend,kbeg,iNOff,ijNOff)
            ijkNode(4) = IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)
          ENDIF  ! lbound 
          DO ic = 1,grid%nCorns(iReg)
            IF (ijkNode(1)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(2)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(3)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(4)==grid%ijkCorn(ic,iReg)) grid%nghbor(3,1:nbour,ic)= 2
          ENDDO
        ENDIF    ! bc_external
      ENDDO      ! iPatch

    ENDIF  ! myProcid
  ENDDO    ! iReg

! third step search for external corners

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(iLev)%patches(iPatch)
        lbound =  patch%lbound
        bcType =  patch%bcType

        CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                        ibeg,iend,jbeg,jend,kbeg,kend )

        IF (patch%bcMotion == BC_EXTERNAL) THEN
          IF (lbound==1 .OR. lbound==2) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)
            ijkNode(3) = IndIJK(ibeg,jend,kend,iNOff,ijNOff)
            ijkNode(4) = IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)
          ELSEIF (lbound==3 .OR. lbound==4) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jbeg,kend,iNOff,ijNOff)
            ijkNode(3) = IndIJK(iend,jbeg,kend,iNOff,ijNOff)
            ijkNode(4) = IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)
          ELSEIF (lbound==5 .OR. lbound==6) THEN
            ijkNode(1) = IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff)
            ijkNode(2) = IndIJK(ibeg,jend,kbeg,iNOff,ijNOff)
            ijkNode(3) = IndIJK(iend,jend,kbeg,iNOff,ijNOff)
            ijkNode(4) = IndIJK(iend,jbeg,kbeg,iNOff,ijNOff)
          ENDIF  ! lbound 
          DO ic = 1,grid%nCorns(iReg)
            IF (ijkNode(1)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(2)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(3)==grid%ijkCorn(ic,iReg) .OR. &
                ijkNode(4)==grid%ijkCorn(ic,iReg)) grid%nghbor(3,1:nbour,ic)= 0
          ENDDO
        ENDIF    ! bc_external
      ENDDO      ! iPatch

!      DO ic = 1,grid%nCorns(iReg)
!        DO k=1,nbour
! ------- inhibitor
!          IF (iReg==70) &
! ------- titan4
!          IF (iReg==120) &
!             write(*,*)iReg,ic,k,edgelen,grid%nghbor(1:3,k,ic)
!        ENDDO
!      ENDDO

    ENDIF  ! myProcid
  ENDDO    ! iReg

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameSrchNeighbors


!******************************************************************************
!
! Purpose: correct closest neighbors
!
! Description: none.
!
! Input: regions = data of current region.
!
! Output: grid%nghbor = corrected neighbouring points identified
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameCorrectNeighbors( regions )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetNodeOffset

  IMPLICIT NONE
#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: i, j, k, ic, iReg, nc, nReg, lc, lReg

! ... local variables
  INTEGER :: iLev, bcType, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ijkNode(4), iNOff, ijNOff, errFl
  REAL(RFREAL) :: edgeLen, ds, du2, duMax
  REAL(RFREAL), POINTER :: xyz(:,:)
  REAL(RFREAL), ALLOCATABLE :: dist(:,:)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameCorrectNeighbors',&
  'RFLO_ModMoveGridNconform4.F90' )

! search for six closest neighbours -------------------------------------------

  ALLOCATE( dist(global%moveGridRegNc,global%nRegions), stat=errFl )
  IF (errFl>0) GOTO 88

  iLev = 1

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

      xyz => regions(iReg)%levels(iLev)%gridOld%xyz

! --- calculate the shortest cell edge

      edgeLen = 1.E+30_RFREAL

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ijkNode(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
            ijkNode(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
            ijkNode(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
            ijkNode(4) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
            ds      = SQRT((xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
            ds      = SQRT((xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
            ds      = SQRT((xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(1)))**2+ &
                           (xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(1)))**2+ &
                           (xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(1)))**2)
            edgeLen = MIN(edgeLen,ds)
          ENDDO
        ENDDO
      ENDDO

      DO ic = 1,grid%nCorns(iReg)
        DO k = 1,global%moveGridNbour
          nc   = grid%nghbor(1,k,ic)
          nReg = grid%nghbor(2,k,ic)
          duMax = -1.e+20_RFREAL

          DO lReg = 1,global%nRegions
            DO lc = 1,grid%nCorns(lReg)
              dist(lc,lReg) = SQRT((grid%regCornOrig(XCOORD,nc,nReg)- &
                                    grid%regCornOrig(XCOORD,lc,lReg))**2 + &
                                    (grid%regCornOrig(YCOORD,nc,nReg)- &
                                    grid%regCornOrig(YCOORD,lc,lReg))**2 + &
                                    (grid%regCornOrig(ZCOORD,nc,nReg)- &
                                    grid%regCornOrig(ZCOORD,lc,lReg))**2)

              IF (dist(lc,lReg) < 0.1_RFREAL*edgelen) THEN
                du2 = grid%regCornOld(XCOORD,lc,lReg)**2 + &
                      grid%regCornOld(YCOORD,lc,lReg)**2 + &
                      grid%regCornOld(ZCOORD,lc,lReg)**2

                IF ( du2 > duMax ) THEN
                  dumax = du2
                  grid%nghbor(1,k,ic) = lc
                  grid%nghbor(2,k,ic) = lReg
                ENDIF  ! duMax
              ENDIF    ! dist
            ENDDO  ! lc
          ENDDO    ! lReg
        ENDDO    ! k
      ENDDO      ! ic

    ENDIF  ! myProcid
  ENDDO    ! iReg

! deallocate temporary arrays -------------------------------------------------

  DEALLOCATE( dist, stat=errFl ); IF (errFl>0) GOTO 99

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameCorrectNeighbors

!******************************************************************************
!
! Purpose: move block corners by averaging over closest neighbours
!
! Description: none.
!
! Input: regions = data of current region.
!
! Output: region%levels%grid%regCorn = new block corners movement.
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameMoveCorners( regions )

  USE ModTools, ONLY : IsNan
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, ico, k, l

! ... local variables
  INTEGER :: iLev, interior, ijkCorn, nbour, errFl
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER, ALLOCATABLE :: nco(:), nReg(:)
  REAL(RFREAL) :: rdenom, amp(3), pow, sum
  REAL(RFREAL), ALLOCATABLE :: dist(:), wght(:)

  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameMoveCorners',&
  'RFLO_ModMoveGridNconform4.F90' )

! move block corners ----------------------------------------------------------

  iLev   = 1
  amp(1) = global%moveGridAmplifX
  amp(2) = global%moveGridAmplifY
  amp(3) = global%moveGridAmplifZ
  pow    = global%moveGridPower
  nbour  = global%moveGridNbour

  ALLOCATE( nco(nbour), stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( nreg(nbour), stat=errFl ); IF (errFl>0) GOTO 88

  ALLOCATE( dist(nbour), stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( wght(nbour), stat=errFl ); IF (errFl>0) GOTO 88

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      DO ico = 1,grid%nCorns(iReg)
        nco(1:nbour)  = grid%nghbor(1,1:nbour,ico)
        nReg(1:nbour) = grid%nghbor(2,1:nbour,ico)
        interior  = grid%nghbor(3,1  ,ico)

        IF (interior==1) THEN
          DO k = 1,nbour
            dist(k) = (grid%regCornOrig(XCOORD,nco(k),nReg(k)) - &
                       grid%regCornOrig(XCOORD,ico,iReg))**2 +   &
                      (grid%regCornOrig(YCOORD,nco(k),nReg(k)) - &
                       grid%regCornOrig(YCOORD,ico,iReg))**2 +   &
                      (grid%regCornOrig(ZCOORD,nco(k),nReg(k)) - &
                       grid%regCornOrig(ZCOORD,ico,iReg))**2
            dist(k) = 1._RFREAL/SQRT( dist(k) )**pow
          ENDDO

          sum = 0._RFREAL
          DO l = 1,nbour
            sum = sum + dist(l)
          ENDDO
          rdenom = 1._RFREAL/sum

          DO k = 1,nbour
            wght(k) = dist(k)*rdenom
!            write(*,*)iReg,ico,k,nReg(k),nco(k),dist(k),wght(k)
            IF (IsNan(wght(k))) &
              CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
                              'invalid weights for global frame motion')
          ENDDO

          sum = 0._RFREAL
          DO l = 1,nbour
            sum = sum + wght(l)*grid%regCornOld(XCOORD,nco(l),nReg(l))
          ENDDO
          grid%regCorn(XCOORD,ico,iReg) = sum

          sum = 0._RFREAL
          DO l = 1,nbour
            sum = sum + wght(l)*grid%regCornOld(YCOORD,nco(l),nReg(l))
          ENDDO
          grid%regCorn(YCOORD,ico,iReg) = sum

          sum = 0._RFREAL
          DO l = 1,nbour
            sum = sum + wght(l)*grid%regCornOld(ZCOORD,nco(l),nReg(l))
          ENDDO
          grid%regCorn(ZCOORD,ico,iReg) = sum

        ENDIF  ! interior
      ENDDO    ! ico
    ENDIF      ! myProcid
  ENDDO        ! iReg

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid => regions(iReg)%levels(iLev)%grid

      DO ico = 1,grid%nCorns(iReg)
        interior = grid%nghbor(3, 1, ico)
        IF (interior==1) THEN
          grid%regCornOld(XCOORD,ico,iReg)=amp(1)*grid%regCorn(XCOORD,ico,iReg)
          grid%regCornOld(YCOORD,ico,iReg)=amp(2)*grid%regCorn(YCOORD,ico,iReg)
          grid%regCornOld(ZCOORD,ico,iReg)=amp(3)*grid%regCorn(ZCOORD,ico,iReg)

          ijkCorn = grid%ijkCorn(ico,iReg)
          grid%xyz(XCOORD,ijkCorn) = grid%regCornOld(XCOORD,ico,iReg)
          grid%xyz(YCOORD,ijkCorn) = grid%regCornOld(YCOORD,ico,iReg)
          grid%xyz(ZCOORD,ijkCorn) = grid%regCornOld(ZCOORD,ico,iReg)
        ENDIF
      ENDDO    ! ico

    ENDIF      ! myProcid
  ENDDO        ! iReg

! deallocate temporary arrays

  DEALLOCATE( nco, stat=errFl );  IF (errFl>0) GOTO 99
  DEALLOCATE( nreg, stat=errFl ); IF (errFl>0) GOTO 99

  DEALLOCATE( dist, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( wght, stat=errFl ); IF (errFl>0) GOTO 99

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameMoveCorners


!******************************************************************************
!
! Purpose: shift corners orthogonally to solid surfaces and averaged over
!          shared neighbors
!
! Description: none.
!
! Input: regions = data of current region.
!
! Output: region%levels%grid%regCorn(Old) = new block corners movement.
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameOrthoShift( regions )

  USE RFLO_ModVectorTensor, ONLY : RFLO_NormCrossProd

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, ico, lb, mc, k

! ... local variables
  INTEGER :: iLev, interior, ijkCorn, ind(6,4,4), kf, nc, nReg
  INTEGER :: lbb, lbe, orthDir
  REAL(RFREAL) :: rlen, eps, orthWg(XCOORD:ZCOORD)
  REAL(RFREAL) :: dif(XCOORD:ZCOORD), shift(XCOORD:ZCOORD)
  REAL(RFREAL) :: s1(XCOORD:ZCOORD), s2(XCOORD:ZCOORD), s3(XCOORD:ZCOORD) 
  REAL(RFREAL), POINTER :: cBuff(:,:,:), cOrth(:,:,:)
  LOGICAL :: solidReg

  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameOrthoShift',&
  'RFLO_ModMoveGridNconform4.F90' )

! set edge ind ----------------------------------------------------------------
! ind(:,:,1)  = opposite corner of current patch corner
! ind(:,:,2)  = current patch corner
! ind(:,:,3)  = left  corner of current patch corner for inward screw
! ind(:,:,4)  = right corner of current patch corner for inward screw

  ind(1, 1,:) = (/5, 1, 4, 2/)
  ind(1, 2,:) = (/6, 2, 1, 3/)
  ind(1, 3,:) = (/7, 3, 2, 4/)
  ind(1, 4,:) = (/8, 4, 3, 1/)
  ind(2, 1,:) = (/1, 5, 6, 8/)
  ind(2, 2,:) = (/2, 6, 7, 5/)
  ind(2, 3,:) = (/3, 7, 8, 6/)
  ind(2, 4,:) = (/4, 8, 5, 7/)
  ind(3, 1,:) = (/4, 1, 2, 5/)
  ind(3, 2,:) = (/8, 5, 1, 6/)
  ind(3, 3,:) = (/7, 6, 5, 2/)
  ind(3, 4,:) = (/3, 2, 6, 1/)
  ind(4, 1,:) = (/1, 4, 8, 3/)
  ind(4, 2,:) = (/5, 8, 7, 4/)
  ind(4, 3,:) = (/6, 7, 3, 8/)
  ind(4, 4,:) = (/2, 3, 4, 7/)
  ind(5, 1,:) = (/2, 1, 5, 4/)
  ind(5, 2,:) = (/3, 4, 1, 8/)
  ind(5, 3,:) = (/7, 8, 4, 5/)
  ind(5, 4,:) = (/6, 5, 8, 1/)
  ind(6, 1,:) = (/1, 2, 3, 6/)
  ind(6, 2,:) = (/4, 3, 7, 2/)
  ind(6, 3,:) = (/8, 7, 6, 3/)
  ind(6, 4,:) = (/5, 6, 2, 7/)

! get parameters --------------------------------------------------------------

  orthDir        = global%moveGridOrthDir
  orthWg(XCOORD) = global%moveGridOrthWghtX *(1._RFREAL + global%skewness)
  orthWg(YCOORD) = global%moveGridOrthWghtY *(1._RFREAL + global%skewness)
  orthWg(ZCOORD) = global%moveGridOrthWghtZ *(1._RFREAL + global%skewness)
                                            ! this to adjust wght by skewness

! move block corners ----------------------------------------------------------

  iLev = 1
  s3   = 0._RFREAL

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid   => regions(iReg)%levels(iLev)%grid
      cBuff  => grid%regCornBuff
      cOrth  => grid%regCornOrth

      cOrth  = cBuff

      solidReg = .FALSE. 
      DO lb = 1,6
        IF (grid%boundMoved(lb)) solidReg = .TRUE.
      ENDDO

      IF (orthDir==OFF) THEN          ! apply orthogonality to all solid surf.
        lbb=1
        lbe=6
      ELSEIF (orthDir==ICOORD) THEN   ! to i-solid surface only
        lbb=1
        lbe=2
      ELSEIF (orthDir==JCOORD) THEN   ! to j-solid surface only
        lbb=3
        lbe=4
      ELSEIF (orthDir==KCOORD) THEN   ! to k-solid surface only
        lbb=5
        lbe=6
      ENDIF

!      IF (solidReg) THEN

      DO lb = lbb,lbe

        IF (grid%boundMoved(lb)) THEN
          DO mc = 1,4
! --------- shift corner-ind(lb,mc,1)

            rlen = SQRT( (cBuff(XCOORD,ind(lb,mc,1),iReg)- &
                          cBuff(XCOORD,ind(lb,mc,2),iReg))**2 + &
                         (cBuff(YCOORD,ind(lb,mc,1),iReg)- &
                          cBuff(YCOORD,ind(lb,mc,2),iReg))**2 + &
                         (cBuff(ZCOORD,ind(lb,mc,1),iReg)- &
                          cBuff(ZCOORD,ind(lb,mc,2),iReg))**2 )
            s1(:) = cBuff(:,ind(lb,mc,3),iReg)-cBuff(:,ind(lb,mc,2),iReg)
            s2(:) = cBuff(:,ind(lb,mc,4),iReg)-cBuff(:,ind(lb,mc,2),iReg)
            CALL RFLO_NormCrossProd( s1,s2,s3 )
            cOrth(:,ind(lb,mc,1),iReg) = cBuff(:,ind(lb,mc,2),iReg)+ &
                                         rlen*s3(:)
          ENDDO ! mc
        ENDIF   ! boundMoved
      ENDDO     ! lb

!      ENDIF     ! solidReg

    ENDIF     ! myProcid
  ENDDO       ! iReg

  CALL RFLO_MgFrameBroadCast( regions,3,1 )

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &    ! region active and
        regions(iReg)%active==ACTIVE) THEN               ! on my processor

      grid   => regions(iReg)%levels(iLev)%grid
      cBuff  => grid%regCornBuff
      cOrth  => grid%regCornOrth

      DO ico = 1,grid%nCorns(iReg)
        interior = grid%nghbor(3, 1, ico)

        IF (interior==1) THEN
          ijkCorn = grid%ijkCorn(ico,iReg)

          shift = 0._RFREAL
          eps = 10._RFREAL*EPSILON( 1._RFREAL )
          kf = 0 
          DO k = 1,grid%nshared(ico)
            nc   = grid%cshared(1,k,ico)
            nReg = grid%cshared(2,k,ico)
            dif(:)   = cOrth(:,nc,nReg) - cBuff(:,nc,nReg)
            IF ((ABS(dif(XCOORD)) > eps) .OR. &
                (ABS(dif(YCOORD)) > eps) .OR. &
                (ABS(dif(ZCOORD)) > eps)) kf = kf+1
            shift(:) = shift(:) + dif(:)
          ENDDO
          IF ((ABS(shift(XCOORD)) > eps) .OR. &
              (ABS(shift(YCOORD)) > eps) .OR. &
              (ABS(shift(ZCOORD)) > eps)) THEN
!            grid%regCornOld(:,ico,iReg) = grid%regCornOld(:,ico,iReg) + &
!                                          orthWg(:)*shift(:)/grid%nshared(ico)
            IF ((ABS(grid%regCornOld(XCOORD,ico,iReg)) > eps) .OR. &
                (ABS(grid%regCornOld(YCOORD,ico,iReg)) > eps) .OR. &
                (ABS(grid%regCornOld(ZCOORD,ico,iReg)) > eps)) THEN
              grid%regCornOld(:,ico,iReg) = (1._RFREAL-orthWg(:))* &
                                            grid%regCornOld(:,ico,iReg) + &
                                            orthWg(:)*shift(:)/kf
            ELSE
              grid%regCornOld(:,ico,iReg) = shift(:)/kf
            ENDIF
          ENDIF  ! ABS(shift)

          ijkCorn = grid%ijkCorn(ico,iReg)
          grid%xyz(XCOORD,ijkCorn) = grid%regCornOld(XCOORD,ico,iReg)
          grid%xyz(YCOORD,ijkCorn) = grid%regCornOld(YCOORD,ico,iReg)
          grid%xyz(ZCOORD,ijkCorn) = grid%regCornOld(ZCOORD,ico,iReg)

        ENDIF
      ENDDO    ! ico

    ENDIF      ! myProcid
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameOrthoShift


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

SUBROUTINE RFLO_MgFrameSurfaces( regions,someMoved,iType )

  USE ModInterfaces, ONLY : RFLO_GetDeformation, RFLO_ArcLengthBounds, &
                            RFLO_EdgeDeformation, RFLO_BoundaryDeformation, &
                            RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes
  IMPLICIT NONE

! ... parameters
  LOGICAL :: someMoved
  INTEGER :: iType

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iter, iPatch, i, j, k, ijkN

! ... local variables
  INTEGER :: iLev, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff
  TYPE(t_grid), POINTER   :: grid, gridOld
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameSurfaces',&
  'RFLO_ModMoveGridNconform4.F90' )

! move grid separately for each region ----------------------------------------

  someMoved = .false.
  iLev  = 1

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE .AND. &            ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN          ! and moving

      grid      => regions(iReg)%levels(iLev)%grid
      gridOld   => regions(iReg)%levels(iLev)%gridOld
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

! broadcast and compute block corners deformation

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

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE .AND. &            ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN          ! and moving

      grid      => regions(iReg)%levels(iLev)%grid
      gridOld   => regions(iReg)%levels(iLev)%gridOld

! --- calculate deformations at remaining edges

!      CALL RFLO_MgFrameEdgesO( regions(iReg),iType,grid%boundMoved, &
!                              grid%allExternal,grid%edgeMoved,grid%arcLen12, &
!                              grid%arcLen34,grid%arcLen56,gridOld%xyzOld,grid%xyz )
      CALL RFLO_MgFrameEdges( regions(iReg),grid%edgeMoved, &
                              gridOld%xyzOld,grid%xyz )

!linMotion      CALL RFLO_MgFrameRestoreExternal( regions(iReg) )

! --- calculate deformations at remaining boundaries

      IF (iType==1) THEN
!        CALL RFLO_MgFrameBndDeformation0( regions(iReg),grid%boundMoved, &
!                                       grid%edgeMoved,grid%arcLen12, &
!                                       grid%arcLen34,grid%arcLen56, &
!                                       gridOld%xyzOld,grid%xyz )
        CALL RFLO_MgFrameBndDeformation( regions(iReg),gridOld%xyzOld,grid%xyz )
                                         
      ELSE  ! iType
        CALL RFLO_BoundaryDeformation( regions(iReg),grid%boundMoved, &
                                       grid%edgeMoved,grid%arcLen12, &
                                       grid%arcLen34,grid%arcLen56, &
                                       gridOld%xyzOld,grid%xyz )
      ENDIF ! iType
!      CALL RFLO_MgFrameRestoreExternal( regions(iReg) )

    ENDIF   ! region on this processor and active, grid moving
  ENDDO     ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameSurfaces


!******************************************************************************
!
! Purpose: restore deformation of solid surfaces from genx at given patches.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%grid%xyz = deformations at the boundaries restored
!
! Notes: grid%xyz temporarily stores nodal displacements. The 'untouched'
!        deformation from genx has been saved in grid%xyzOld.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameRestoreExternal( region )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes
  IMPLICIT NONE
#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iReg, iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, ijkN, lbound
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff
  TYPE(t_grid), POINTER   :: grid
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_MgFrameRestoreExternal',&
  'RFLO_ModMoveGridNconform4.F90' )

! parameters and pointers -----------------------------------------------------

  iLev = 1
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  grid => region%levels(iLev)%grid

! restore displacements

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    lbound =  patch%lbound

    IF (patch%bcMotion == BC_EXTERNAL .AND. &
        (grid%allExternal(lbound).EQV..FALSE.)) THEN
!    IF (patch%bcMotion == BC_EXTERNAL) THEN
      CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i,j,k,iNOff,ijNOff)
            grid%xyz(XCOORD,ijkN) = grid%xyzOld(XCOORD,ijkN)
            grid%xyz(YCOORD,ijkN) = grid%xyzOld(YCOORD,ijkN)
            grid%xyz(ZCOORD,ijkN) = grid%xyzOld(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDDO

    ENDIF    ! external BC
  ENDDO      ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MgFrameRestoreExternal

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

SUBROUTINE RFLO_MgFrameEdgesO( region,iType,boundMoved,allExternal,edgeMoved, &
                              arcLen12,arcLen34,arcLen56,xyzOld,dNode )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint1d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), allExternal(6), edgeMoved(12)

  INTEGER :: iType
  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)

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

  CALL RegisterFunction( region%global,'RFLO_MgFrameEdgesO',&
  'RFLO_ModMoveGridNconform4.F90' )

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

  edgeMoved(:) = .false.

  IF (iType/=1) THEN
    IF (boundMoved(1) .AND. allExternal(1)) THEN
      edgeMoved( 1) = .true.; edgeMoved( 2) = .true.
      edgeMoved( 3) = .true.; edgeMoved( 4) = .true.
    ENDIF
    IF (boundMoved(2) .AND. allExternal(2)) THEN
      edgeMoved( 5) = .true.; edgeMoved( 6) = .true.
      edgeMoved( 7) = .true.; edgeMoved( 8) = .true.
    ENDIF
    IF (boundMoved(3) .AND. allExternal(3)) THEN
      edgeMoved( 1) = .true.; edgeMoved( 5) = .true.
      edgeMoved( 9) = .true.; edgeMoved(10) = .true.
    ENDIF
    IF (boundMoved(4) .AND. allExternal(4)) THEN
      edgeMoved( 3) = .true.; edgeMoved( 7) = .true.
      edgeMoved(11) = .true.; edgeMoved(12) = .true.
    ENDIF
    IF (boundMoved(5) .AND. allExternal(5)) THEN
      edgeMoved( 4) = .true.; edgeMoved( 8) = .true.
      edgeMoved( 9) = .true.; edgeMoved(11) = .true.
    ENDIF
    IF (boundMoved(6) .AND. allExternal(6)) THEN
      edgeMoved( 2) = .true.; edgeMoved( 6) = .true.
      edgeMoved(10) = .true.; edgeMoved(12) = .true.
    ENDIF
  ENDIF  ! iType

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
            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
          ELSE IF (switch(iEdge,5) == 34) THEN
            ijkN     = IndIJK(l2c,ind   ,l1c,iNOff,ijNOff)
            ijkN1    = IndIJK(l2c,ind-1 ,l1c,iNOff,ijNOff)
            ijkNBeg  = IndIJK(l2c,indBeg,l1c,iNOff,ijNOff)
            ijkNEnd  = IndIJK(l2c,indEnd,l1c,iNOff,ijNOff)
            arcLen   = arcLen34(l1c,l2c)
            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
          ELSE IF (switch(iEdge,5) == 56) THEN
            ijkN     = IndIJK(l1c,l2c,ind   ,iNOff,ijNOff)
            ijkN1    = IndIJK(l1c,l2c,ind-1 ,iNOff,ijNOff)
            ijkNBeg  = IndIJK(l1c,l2c,indBeg,iNOff,ijNOff)
            ijkNEnd  = IndIJK(l1c,l2c,indEnd,iNOff,ijNOff)
            arcLen   = arcLen56(l1c,l2c)
            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
          ENDIF
          ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                         (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                         (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
          s  = ds/arcLen

          CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
          dNode(:,ijkN) = dN(:)
        ENDDO   ! i
      ENDIF   ! nghbor
    ENDIF     ! edgeMoved
  ENDDO       ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgFrameEdgesO

!******************************************************************************
!
! Purpose: calculate node displacements on those patch edges whose end points 
!          both are or one of them is interior point.
!
! Description: points along patch edge are shifted using 1-D linear transfinite
!              interpolation (TFI).
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        xyzOld     = grid from previous time step.
!
! Output: edgeMoved = flag if discretization at an edge was changed
!         dNode     = updated deformations at edges.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameEdges( region,edgeMoved,xyzOld,dNode )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, &
                            RFLO_Tfint1d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  LOGICAL :: edgeMoved(12)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)

! ... loop variables
  INTEGER :: ipEdge, ind, iPatch, ic

! ... local variables
  INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, ib, ie, jb, je , kb, ke
  INTEGER :: indBeg, indEnd, ijkN, ijkN1, ijkNBeg, ijkNEnd, iNOff, ijNOff
  INTEGER :: lbound, intb, inte, iReg, iEdge, interType, ijkTest

  REAL(RFREAL)  :: arcLen, ds, s, dN(3), dNBeg(3), dNEnd(3)
  TYPE(t_patch), POINTER :: patch
  TYPE(t_grid), POINTER  :: grid

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MgFrameEdges',&
  'RFLO_ModMoveGridNconform4.F90' )

! get dimensions and pointers -------------------------------------------------

  iReg = region%iRegionGlobal
  iLev = 1
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  grid => region%levels(iLev)%grid

! edge movement flag ----------------------------------------------------------

  edgeMoved(:) = .true.

! loop over all patch edges ---------------------------------------------------

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    lbound =  patch%lbound

    CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                    ibeg,iend,jbeg,jend,kbeg,kend )


!linMotion    IF (patch%bcMotion == BC_EXTERNAL) GOTO 777
    DO ipEdge = 1,4

      IF (lbound==1 .OR. lbound==2) THEN
        IF (ipEdge==1) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jbeg
          je       = jbeg
          kb       = kbeg
          ke       = kend
          indBeg   = kb
          indEnd   = ke
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==1) THEN
            iEdge = 1
          ELSE
            iEdge = 5
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(1) = intb
          patch%position(2) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN
               
            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,jbeg,ind  ,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,jbeg,ind-1,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,jbeg,ind  ,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,jbeg,ind-1,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==2) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jbeg
          je       = jend
          kb       = kend
          ke       = kend
          indBeg   = jb
          indEnd   = je
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==1) THEN
            iEdge = 2
          ELSE
            iEdge = 6
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
!            IF (region%iRegionGlobal==70) write(*,*)iReg,ic, &
!              grid%ijkCorn(ic,iReg),ijkNBeg,ijkNEnd,intb,inte, &
!              grid%nghbor(3,1,ic)
          ENDDO
          patch%position(2) = intb
          patch%position(3) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
      
            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,ind  ,kend,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,ind-1,kend,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,ind  ,kend,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,ind-1,kend,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==3) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jend
          je       = jend
          kb       = kbeg
          ke       = kend
          indBeg   = kb
          indEnd   = ke
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==1) THEN
            iEdge = 3
          ELSE
            iEdge = 7
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(3) = inte
          patch%position(4) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
      
            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,jend,ind  ,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,jend,ind-1,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,jend,ind  ,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,jend,ind-1,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==4) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jbeg
          je       = jend
          kb       = kbeg
          ke       = kbeg
          indBeg   = jb
          indEnd   = je
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==1) THEN
            iEdge = 4
          ELSE
            iEdge = 8
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(4) = inte
          patch%position(1) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
      
            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,ind  ,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,ind-1,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,ind  ,kbeg,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,ind-1,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
!              IF (region%iRegionGlobal==10 .AND. lbound==2) &
!                 write(*,*)ind,s,arcLen,dNBeg,dNEnd
            ENDDO
!linMotion          ENDIF  ! intb
        ENDIF  ! ipEdge
      ENDIF  ! lbound

      IF (lbound==3 .OR. lbound==4) THEN
        IF (ipEdge==1) THEN
          ib       = ibeg
          ie       = iend
          jb       = jbeg
          je       = jbeg
          kb       = kbeg
          ke       = kbeg
          indBeg   = ib
          indEnd   = ie
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==3) THEN
            iEdge = 9
          ELSE
            iEdge = 12
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(1) = intb
          patch%position(2) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ind  ,jbeg,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jbeg,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ind  ,jbeg,kbeg,iNOff,ijNOff)
              ijkN1    = IndIJK(ind-1,jbeg,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==2) THEN
          ib       = iend
          ie       = iend
          jb       = jbeg
          je       = jbeg
          kb       = kbeg
          ke       = kend
          indBeg   = kb
          indEnd   = ke
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==3) THEN
            iEdge = 5
          ELSE
            iEdge = 7
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(2) = intb
          patch%position(3) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
      
            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(iend,jbeg,ind  ,iNOff,ijNOff)
              ijkN1  = IndIJK(iend,jbeg,ind-1,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(iend,jbeg,ind  ,iNOff,ijNOff)
              ijkN1    = IndIJK(iend,jbeg,ind-1,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==3) THEN
          ib       = ibeg
          ie       = iend
          jb       = jbeg
          je       = jbeg
          kb       = kend
          ke       = kend
          indBeg   = ib
          indEnd   = ie
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==3) THEN
            iEdge = 10
          ELSE
            iEdge = 11
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(3) = inte
          patch%position(4) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion             ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ind  ,jbeg,kend,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jbeg,kend,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ind  ,jbeg,kend,iNOff,ijNOff)
              ijkN1    = IndIJK(ind-1,jbeg,kend,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==4) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jbeg
          je       = jbeg
          kb       = kbeg
          ke       = kend
          indBeg   = kb
          indEnd   = ke
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==3) THEN
            iEdge = 1
          ELSE
            iEdge = 3
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(4) = inte
          patch%position(1) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)
      
            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,jbeg,ind  ,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,jbeg,ind-1,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,jbeg,ind  ,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,jbeg,ind-1,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF  ! intb
        ENDIF  ! ipEdge
      ENDIF  ! lbound

      IF (lbound==5 .OR. lbound==6) THEN
        IF (ipEdge==1) THEN
          ib       = ibeg
          ie       = ibeg
          jb       = jbeg
          je       = jend
          kb       = kbeg
          ke       = kbeg
          indBeg   = jb
          indEnd   = je
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==5) THEN
            iEdge = 4
          ELSE
            iEdge = 2
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(1) = intb
          patch%position(2) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ibeg,ind  ,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ibeg,ind-1,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(ibeg,ind  ,kbeg,iNOff,ijNOff)
              ijkN1    = IndIJK(ibeg,ind-1,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==2) THEN
          ib       = ibeg
          ie       = iend
          jb       = jend
          je       = jend
          kb       = kbeg
          ke       = kbeg
          indBeg   = ib
          indEnd   = ie
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==5) THEN
            iEdge = 12
          ELSE
            iEdge = 11
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(2) = intb
          patch%position(3) = inte

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ind  ,jend,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jend,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN   = IndIJK(ind  ,jend,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jend,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==3) THEN
          ib       = iend
          ie       = iend
          jb       = jbeg
          je       = jend
          kb       = kbeg
          ke       = kbeg
          indBeg   = jb
          indEnd   = je
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==5) THEN
            iEdge = 8
          ELSE
            iEdge = 6
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(3) = inte
          patch%position(4) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(iend,ind  ,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(iend,ind-1,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN     = IndIJK(iend,ind  ,kbeg,iNOff,ijNOff)
              ijkN1    = IndIJK(iend,ind-1,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF
        ELSEIF (ipEdge==4) THEN
          ib       = ibeg
          ie       = iend
          jb       = jbeg
          je       = jbeg
          kb       = kbeg
          ke       = kbeg
          indBeg   = ib
          indEnd   = ie
          ijkNBeg  = IndIJK(ib,jb,kb,iNOff,ijNOff)
          ijkNEnd  = IndIJK(ie,je,ke,iNOff,ijNOff)
          ijkTest  = IndIJK((ib+ie)/2,(jb+je)/2,(kb+ke)/2,iNOff,ijNOff)
          intb     = 1
          inte     = 1
          IF (lbound==5) THEN
            iEdge = 9
          ELSE
            iEdge = 10
          ENDIF
          interType =  region%levels(iLev)%edgeCells(iEdge)%interType
          DO ic = 1,grid%nCorns(iReg)
            IF (grid%ijkCorn(ic,iReg)==ijkNBeg) intb = grid%nghbor(3,1,ic) 
            IF (grid%ijkCorn(ic,iReg)==ijkNEnd) inte = grid%nghbor(3,1,ic)
          ENDDO
          patch%position(4) = inte
          patch%position(1) = intb

!linMotion          IF (((intb/=0 .OR. inte/=0) .AND. (interType==EDGE_INTERACT_FULL .AND. &
!linMotion                region%levels(iLev)%edgeCells(iEdge)%interact .EQV. .true.)) &
!linMotion                         .OR. &
!linMotion!              (intb==2 .OR. inte==2)) THEN
!linMotion              (intb==2 .OR. inte==2) &
!linMotion                         .OR. &
!linMotion              ((intb==0 .AND. inte==0) .AND. ABS(dNode(XCOORD,ijkTest)+ &
!linMotion                dNode(YCOORD,ijkTest)+dNode(ZCOORD,ijkTest)) < &
!linMotion                100._RFREAL*EPSILON(1._RFREAL))) THEN

            dNBeg(:) = dNode(:,ijkNBeg)
            dNEnd(:) = dNode(:,ijkNEnd)

            arcLen = 0._RFREAL
            DO ind=indBeg+1,indEnd
              ijkN   = IndIJK(ind  ,jbeg,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jbeg,kbeg,iNOff,ijNOff)
              arcLen = arcLen + &
                       SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                            (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                            (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
            ENDDO
            ds = 0._RFREAL
            DO ind=indBeg+1,indEnd-1
              ijkN   = IndIJK(ind  ,jbeg,kbeg,iNOff,ijNOff)
              ijkN1  = IndIJK(ind-1,jbeg,kbeg,iNOff,ijNOff)
              ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                             (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                             (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
              s  = ds/arcLen
              CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
              dNode(:,ijkN) = dN(:)
            ENDDO
!linMotion          ENDIF  ! intb
        ENDIF  ! ipEdge
      ENDIF  ! lbound

    ENDDO    ! ipEdge

777 CONTINUE

  ENDDO      ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgFrameEdges

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

SUBROUTINE RFLO_MgFrameInterfaces( regions,iType )

  USE ModInterfaces, ONLY : RFLO_ExchangeDnodeCopy, RFLO_EdgeDeformation, &
        RFLO_BoundaryDeformation, RFLO_ExchangeDnodeSend, &
        RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER ::iType

! ... loop variables
  INTEGER :: iReg, iPatch, iPass

! ... local variables
  INTEGER :: bcType, iRegSrc, iPatchSrc, lType, nPass

  TYPE(t_grid), POINTER   :: grid, gridOld, gridSrc
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MgFrameInterfaces',&
  'RFLO_ModMoveGridNconform4.F90' )

! fix interfaces between regions ----------------------------------------------

  lType = 1
  nPass = global%moveGridNsmatch
  nPass = MAX( nPass,4 )
  DO iPass=1,nPass

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
              IF (iPass < nPass-1 .AND. lType==1) THEN
                CALL RFLO_MgFrameEdges( regions(iReg),grid%edgeMoved, &
                                        gridOld%xyzOld,grid%xyz )
!linMotion                CALL RFLO_MgFrameRestoreExternal( regions(iReg) )
                CALL RFLO_MgFrameBndDeformation( regions(iReg),gridOld%xyzOld, &
                                        grid%xyz )
              ELSEIF (iPass < nPass-1 .AND. lType==2) THEN
                CALL RFLO_MgFrameEdgesO( regions(iReg),2,grid%boundMoved, &
                                        grid%allExternal,grid%edgeMoved, &
                                        grid%arcLen12,grid%arcLen34, &
                                        grid%arcLen56,gridOld%xyzOld,grid%xyz )
                CALL RFLO_BoundaryDeformation( regions(iReg), &
                                        grid%boundMoved, &
                                        grid%edgeMoved,grid%arcLen12, &
                                        grid%arcLen34,grid%arcLen56, &
                                        gridOld%xyzOld,grid%xyz )
              ENDIF
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

              IF (iPass < nPass-1 .AND. lType==1) THEN
                CALL RFLO_MgFrameEdges( regions(iReg),grid%edgeMoved, &
                                        gridOld%xyzOld,grid%xyz )
!linMotion                CALL RFLO_MgFrameRestoreExternal( regions(iReg) )
                CALL RFLO_MgFrameBndDeformation( regions(iReg),gridOld%xyzOld, &
                                        grid%xyz )
              ELSEIF (iPass < nPass-1 .AND. lType==2) THEN
                CALL RFLO_MgFrameEdgesO( regions(iReg),2,grid%boundMoved, &
                                        grid%allExternal,grid%edgeMoved, &
                                        grid%arcLen12,grid%arcLen34, &
                                        grid%arcLen56,gridOld%xyzOld,grid%xyz )
                CALL RFLO_BoundaryDeformation( regions(iReg), &
                                        grid%boundMoved, &
                                        grid%edgeMoved,grid%arcLen12, &
                                        grid%arcLen34,grid%arcLen56, &
                                        gridOld%xyzOld,grid%xyz )
              ENDIF
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

END SUBROUTINE RFLO_MgFrameInterfaces

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
!        xyzOld     = grid from previous time step.
!
! Output: dNode = updated deformations at boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameBndDeformation0( region,boundMoved,edgeMoved, &
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

  CALL RegisterFunction( region%global,'RFLO_MgFrameBndDeformation0',&
  'RFLO_ModMoveGridNconform4.F90' )

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
!    IF ((.NOT.boundMoved(iBound)) .AND. &
!        (edgeMoved(switch(iBound,1)) .OR. edgeMoved(switch(iBound,2)) .OR. &
!         edgeMoved(switch(iBound,3)) .OR. edgeMoved(switch(iBound,4)))) THEN

    IF ((edgeMoved(switch(iBound,1)) .OR. edgeMoved(switch(iBound,2)) .OR. &
         edgeMoved(switch(iBound,3)) .OR. edgeMoved(switch(iBound,4)))) THEN

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

END SUBROUTINE RFLO_MgFrameBndDeformation0

!******************************************************************************
!
! Purpose: calculate node displacements on non-external patches
!          (finest grid only).
!
! Description: none.
!
! Input: region     = grid dimensions
!        xyzOld     = grid from previous time step.
!
! Output: dNode = updated deformations at boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************

SUBROUTINE RFLO_MgFrameBndDeformation( region,xyzOld,dNode )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, &
                            RFLO_Tfint2d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters

  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, l1, l2

! ... local variables
  INTEGER :: iLev, iBound, ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: l1b, l1e, l2b, l2e, lc, ijkN, ijkE(4), ijkEm(4), iNOff, ijNOff
  INTEGER :: switch(6,9)

  LOGICAL :: sum12
  REAL(RFREAL) :: arcLen(4), ds(4), s(4)
  REAL(RFREAL) :: e1(3), e2(3), e3(3), e4(3), &
                  p1(3), p2(3), p3(3), p4(3), dN(3)
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MgFrameBndDeformation',&
  'RFLO_ModMoveGridNconform4.F90' )

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set boundary switch ---------------------------------------------------------
! switch(:,1-4) = numbers of the 4 edges of a boundary
! switch(:,5-6) = first/last index in l1-direction
! switch(:,7-8) = first/last index in l2-direction
! switch(:,  9) = constant index

! move nodes on boundaries with active edges ----------------------------------

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    iBound =  patch%lbound

    CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                    ibeg,iend,jbeg,jend,kbeg,kend )

!linMotion    IF (patch%bcMotion/=BC_EXTERNAL) THEN

      switch(1,:) = (/ 1,  2,  3,  4, jbeg, jend, kbeg, kend, ibeg/)
      switch(2,:) = (/ 5,  6,  7,  8, jbeg, jend, kbeg, kend, iend/)
      switch(3,:) = (/ 1,  5,  9, 10, kbeg, kend, ibeg, iend, jbeg/)
      switch(4,:) = (/ 3,  7, 11, 12, kbeg, kend, ibeg, iend, jend/)
      switch(5,:) = (/ 4,  8,  9, 11, ibeg, iend, jbeg, jend, kbeg/)
      switch(6,:) = (/ 2,  6, 10, 12, ibeg, iend, jbeg, jend, kend/)

      l1b = switch(iBound,5)
      l1e = switch(iBound,6)
      l2b = switch(iBound,7)
      l2e = switch(iBound,8)
      lc  = switch(iBound,9)

      p1(:) = dNode(:,patch%corns(1))
      p2(:) = dNode(:,patch%corns(4))
      p3(:) = dNode(:,patch%corns(3))
      p4(:) = dNode(:,patch%corns(2))

! --- compute arclen along patch edges

      IF (iBound==1 .OR. iBound==2) THEN
        arclen(1:2) = 0._RFREAL
        DO l2=l2b+1,l2e
          ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
          ijkE(1)   = IndIJK(lc,jbeg  ,l2    ,iNOff,ijNOff)
          ijkEm(1)  = IndIJK(lc,jbeg  ,l2-1  ,iNOff,ijNOff)
          ijkE(2)   = IndIJK(lc,jend  ,l2    ,iNOff,ijNOff)
          ijkEm(2)  = IndIJK(lc,jend  ,l2-1  ,iNOff,ijNOff)
          arclen(1) = arclen(1) + &
                SQRT((xyzOld(XCOORD,ijkE(1))-xyzOld(XCOORD,ijkEm(1)))**2 + &
                     (xyzOld(YCOORD,ijkE(1))-xyzOld(YCOORD,ijkEm(1)))**2 + &
                     (xyzOld(ZCOORD,ijkE(1))-xyzOld(ZCOORD,ijkEm(1)))**2)
          arclen(2) = arclen(2) + &
                SQRT((xyzOld(XCOORD,ijkE(2))-xyzOld(XCOORD,ijkEm(2)))**2 + &
                     (xyzOld(YCOORD,ijkE(2))-xyzOld(YCOORD,ijkEm(2)))**2 + &
                     (xyzOld(ZCOORD,ijkE(2))-xyzOld(ZCOORD,ijkEm(2)))**2)
        ENDDO
        arclen(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e
          ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
          ijkE(3)   = IndIJK(lc,l1    ,kbeg  ,iNOff,ijNOff)
          ijkEm(3)  = IndIJK(lc,l1-1  ,kbeg  ,iNOff,ijNOff)
          ijkE(4)   = IndIJK(lc,l1    ,kend  ,iNOff,ijNOff)
          ijkEm(4)  = IndIJK(lc,l1-1  ,kend  ,iNOff,ijNOff)
          arclen(3) = arclen(3) + &
                SQRT((xyzOld(XCOORD,ijkE(3))-xyzOld(XCOORD,ijkEm(3)))**2 + &
                     (xyzOld(YCOORD,ijkE(3))-xyzOld(YCOORD,ijkEm(3)))**2 + &
                     (xyzOld(ZCOORD,ijkE(3))-xyzOld(ZCOORD,ijkEm(3)))**2)
          arclen(4) = arclen(4) + &
                SQRT((xyzOld(XCOORD,ijkE(4))-xyzOld(XCOORD,ijkEm(4)))**2 + &
                     (xyzOld(YCOORD,ijkE(4))-xyzOld(YCOORD,ijkEm(4)))**2 + &
                     (xyzOld(ZCOORD,ijkE(4))-xyzOld(ZCOORD,ijkEm(4)))**2)
        ENDDO
      ELSE IF (iBound==3 .OR. iBound==4) THEN
        arclen(1:2) = 0._RFREAL
        DO l2=l2b+1,l2e
          ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
          ijkE(1)   = IndIJK(l2    ,lc,kbeg  ,iNOff,ijNOff)
          ijkEm(1)  = IndIJK(l2-1  ,lc,kbeg  ,iNOff,ijNOff)
          ijkE(2)   = IndIJK(l2    ,lc,kend  ,iNOff,ijNOff)
          ijkEm(2)  = IndIJK(l2-1  ,lc,kend  ,iNOff,ijNOff)
          arclen(1) = arclen(1) + &
                SQRT((xyzOld(XCOORD,ijkE(1))-xyzOld(XCOORD,ijkEm(1)))**2 + &
                     (xyzOld(YCOORD,ijkE(1))-xyzOld(YCOORD,ijkEm(1)))**2 + &
                     (xyzOld(ZCOORD,ijkE(1))-xyzOld(ZCOORD,ijkEm(1)))**2)
          arclen(2) = arclen(2) + &
                SQRT((xyzOld(XCOORD,ijkE(2))-xyzOld(XCOORD,ijkEm(2)))**2 + &
                     (xyzOld(YCOORD,ijkE(2))-xyzOld(YCOORD,ijkEm(2)))**2 + &
                     (xyzOld(ZCOORD,ijkE(2))-xyzOld(ZCOORD,ijkEm(2)))**2)
        ENDDO
        arclen(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e
          ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
          ijkE(3)   = IndIJK(ibeg  ,lc,l1    ,iNOff,ijNOff)
          ijkEm(3)  = IndIJK(ibeg  ,lc,l1-1  ,iNOff,ijNOff)
          ijkE(4)   = IndIJK(iend  ,lc,l1    ,iNOff,ijNOff)
          ijkEm(4)  = IndIJK(iend  ,lc,l1-1  ,iNOff,ijNOff)
          arclen(3) = arclen(3) + &
                SQRT((xyzOld(XCOORD,ijkE(3))-xyzOld(XCOORD,ijkEm(3)))**2 + &
                     (xyzOld(YCOORD,ijkE(3))-xyzOld(YCOORD,ijkEm(3)))**2 + &
                     (xyzOld(ZCOORD,ijkE(3))-xyzOld(ZCOORD,ijkEm(3)))**2)
          arclen(4) = arclen(4) + &
                SQRT((xyzOld(XCOORD,ijkE(4))-xyzOld(XCOORD,ijkEm(4)))**2 + &
                     (xyzOld(YCOORD,ijkE(4))-xyzOld(YCOORD,ijkEm(4)))**2 + &
                     (xyzOld(ZCOORD,ijkE(4))-xyzOld(ZCOORD,ijkEm(4)))**2)
        ENDDO
      ELSE IF (iBound==5 .OR. iBound==6) THEN
        arclen(1:2) = 0._RFREAL
        DO l2=l2b+1,l2e
          ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
          ijkE(1)   = IndIJK(ibeg  ,l2    ,lc,iNOff,ijNOff)
          ijkEm(1)  = IndIJK(ibeg  ,l2-1  ,lc,iNOff,ijNOff)
          ijkE(2)   = IndIJK(iend  ,l2    ,lc,iNOff,ijNOff)
          ijkEm(2)  = IndIJK(iend  ,l2-1  ,lc,iNOff,ijNOff)
          arclen(1) = arclen(1) + &
                SQRT((xyzOld(XCOORD,ijkE(1))-xyzOld(XCOORD,ijkEm(1)))**2 + &
                     (xyzOld(YCOORD,ijkE(1))-xyzOld(YCOORD,ijkEm(1)))**2 + &
                     (xyzOld(ZCOORD,ijkE(1))-xyzOld(ZCOORD,ijkEm(1)))**2)
          arclen(2) = arclen(2) + &
                SQRT((xyzOld(XCOORD,ijkE(2))-xyzOld(XCOORD,ijkEm(2)))**2 + &
                     (xyzOld(YCOORD,ijkE(2))-xyzOld(YCOORD,ijkEm(2)))**2 + &
                     (xyzOld(ZCOORD,ijkE(2))-xyzOld(ZCOORD,ijkEm(2)))**2)
        ENDDO
        arclen(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e
          ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
          ijkE(3)   = IndIJK(l1    ,jbeg  ,lc,iNOff,ijNOff)
          ijkEm(3)  = IndIJK(l1-1  ,jbeg  ,lc,iNOff,ijNOff)
          ijkE(4)   = IndIJK(l1    ,jend  ,lc,iNOff,ijNOff)
          ijkEm(4)  = IndIJK(l1-1  ,jend  ,lc,iNOff,ijNOff)
          arclen(3) = arclen(3) + &
                SQRT((xyzOld(XCOORD,ijkE(3))-xyzOld(XCOORD,ijkEm(3)))**2 + &
                     (xyzOld(YCOORD,ijkE(3))-xyzOld(YCOORD,ijkEm(3)))**2 + &
                     (xyzOld(ZCOORD,ijkE(3))-xyzOld(ZCOORD,ijkEm(3)))**2)
          arclen(4) = arclen(4) + &
                SQRT((xyzOld(XCOORD,ijkE(4))-xyzOld(XCOORD,ijkEm(4)))**2 + &
                     (xyzOld(YCOORD,ijkE(4))-xyzOld(YCOORD,ijkEm(4)))**2 + &
                     (xyzOld(ZCOORD,ijkE(4))-xyzOld(ZCOORD,ijkEm(4)))**2)
        ENDDO
      ENDIF

! --- conduct TFI on interior patch surface

      ds(1:2) = 0._RFREAL
      DO l2=l2b+1,l2e-1

        sum12   = .true.
        ds(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e-1
          IF (iBound==1 .OR. iBound==2) THEN
            ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(lc,jbeg  ,l2    ,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(lc,jbeg  ,l2-1  ,iNOff,ijNOff)
            ijkE(2)   = IndIJK(lc,jend  ,l2    ,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(lc,jend  ,l2-1  ,iNOff,ijNOff)
            ijkE(3)   = IndIJK(lc,l1    ,kbeg  ,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(lc,l1-1  ,kbeg  ,iNOff,ijNOff)
            ijkE(4)   = IndIJK(lc,l1    ,kend  ,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(lc,l1-1  ,kend  ,iNOff,ijNOff)
          ELSE IF (iBound==3 .OR. iBound==4) THEN
            ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(l2    ,lc,kbeg  ,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(l2-1  ,lc,kbeg  ,iNOff,ijNOff)
            ijkE(2)   = IndIJK(l2    ,lc,kend  ,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(l2-1  ,lc,kend  ,iNOff,ijNOff)
            ijkE(3)   = IndIJK(ibeg  ,lc,l1    ,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(ibeg  ,lc,l1-1  ,iNOff,ijNOff)
            ijkE(4)   = IndIJK(iend  ,lc,l1    ,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(iend  ,lc,l1-1  ,iNOff,ijNOff)
          ELSE IF (iBound==5 .OR. iBound==6) THEN
            ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
            ijkE(1)   = IndIJK(ibeg  ,l2    ,lc,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(ibeg  ,l2-1  ,lc,iNOff,ijNOff)
            ijkE(2)   = IndIJK(iend  ,l2    ,lc,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(iend  ,l2-1  ,lc,iNOff,ijNOff)
            ijkE(3)   = IndIJK(l1    ,jbeg  ,lc,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(l1-1  ,jbeg  ,lc,iNOff,ijNOff)
            ijkE(4)   = IndIJK(l1    ,jend  ,lc,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(l1-1  ,jend  ,lc,iNOff,ijNOff)
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

!linMotion    ENDIF      ! not.external
  ENDDO        ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MgFrameBndDeformation

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModMoveGridFrame

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModMoveGridNconform4.F90,v $
! Revision 1.13  2009/08/27 14:04:51  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.12  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/03/18 13:26:00  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.9  2006/03/18 11:03:24  wasistho
! screen printed global skewness and minvol
!
! Revision 1.8  2006/03/16 08:29:56  wasistho
! modified global skewness effect on orthowght
!
! Revision 1.7  2006/03/15 06:38:45  wasistho
! added region and global skewness
!
! Revision 1.6  2006/03/05 22:28:08  wasistho
! fixed syntax error
!
! Revision 1.5  2006/03/05 21:52:28  wasistho
! changed computational space coordinates to be based on initial grid
!
! Revision 1.4  2006/02/11 03:53:30  wasistho
! made some routines public
!
! Revision 1.3  2006/01/28 22:53:00  wasistho
! fixed iEdgeGlo
!
! Revision 1.2  2005/12/05 10:51:19  wasistho
! moved RFLO_NormCrossProd to RFLO_ModVectorTensor
!
! Revision 1.1  2005/09/13 00:07:24  wasistho
! initial import linear approx. solid motion
!
! Revision 1.5  2005/09/09 03:24:58  wasistho
! added lb=1,2 in OrthoShift for more robust option
!
! Revision 1.4  2005/09/05 02:25:22  wasistho
! make 2 passes for ExchangeDnode.. at last
!
! Revision 1.3  2005/08/31 19:01:55  wasistho
! last modified OrthoShift
!
! Revision 1.2  2005/08/28 23:49:30  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.1  2005/08/25 23:11:09  wasistho
! initial import featuring block orthogonality
!
! Revision 1.21  2005/08/18 19:51:27  wasistho
! added user define nPass in mgFrameInterface
!
! Revision 1.20  2005/07/10 21:13:13  wasistho
! global alloc dist in SrchNeighbors and CorrectNeighbors, and added pointer grid in mgFrameMoveCorners
!
! Revision 1.19  2005/07/04 22:07:37  wasistho
! ABS operator was missing in the previous modification
!
! Revision 1.18  2005/07/04 12:07:05  wasistho
! another mgFrameEdge modification to better accomodate non-conforming blocking
!
! Revision 1.17  2005/07/02 23:10:53  wasistho
! modified mgFrameEdges, Surfaces and Interfaces to move edges between slots
!
! Revision 1.16  2005/06/30 19:10:07  wasistho
! made official last added conditions in mgFrameEdges
!
! Revision 1.15  2005/06/30 07:57:00  wasistho
! additional condition to move patchEdges in mgFrameEdges
!
! Revision 1.14  2005/06/29 22:53:31  wasistho
! added interType condition in mgFrameEdges
!
! Revision 1.13  2005/06/29 08:44:03  wasistho
! changed RestoreExternal and set lTyp=1
!
! Revision 1.12  2005/06/29 04:55:48  wasistho
! modified RFLO_MgFrameInterfaces
!
! Revision 1.11  2005/06/27 19:22:50  wasistho
! only ipass=1 in mgFrameInterfaces
!
! Revision 1.10  2005/06/27 01:00:45  wasistho
! stored tolerance in tol
!
! Revision 1.9  2005/06/27 00:36:42  wasistho
! change tolerance in mgFrameSrchNeighbors from 1.e-20 to 1.e-5*edgeLen
!
! Revision 1.8  2005/06/26 06:25:32  wasistho
! nReg==72 to nReg==71
!
! Revision 1.7  2005/06/26 06:11:48  wasistho
! adding more reegions for titan check in mgFrameSrchCorners
!
! Revision 1.6  2005/06/26 05:39:29  wasistho
! bugs fixed in mgFrameCornPoints and mgFrameSrchNeighbors
!
! Revision 1.5  2005/06/25 08:12:19  wasistho
! bug fixed in receiving rvar in mgFrameBCast
!
! Revision 1.4  2005/06/25 06:16:05  wasistho
! swap ENDDO and ENDIF in mgframeBroadCast
!
! Revision 1.3  2005/06/25 03:14:06  wasistho
! enabled nRegions /= nProcs in type 2 gridmotion
!
! Revision 1.2  2005/06/23 08:54:46  wasistho
! fixed bug rdenom should be sum over nbour
!
! Revision 1.1  2005/06/23 01:37:33  wasistho
! make number of closest neighbors user input
!
! Revision 1.12  2005/06/19 12:36:00  wasistho
! update in mgFrameInterfaces
!
! Revision 1.11  2005/06/17 03:51:22  wasistho
! argument in RFLO_exchangeDnodeCopy/Recv is set from true to false
!
! Revision 1.10  2005/06/16 01:28:02  wasistho
! modified mgFrameBndDeformation to allow TFI per patch i.o per block-side
!
! Revision 1.9  2005/06/15 19:20:47  wasistho
! simplified mgFrameInterfaces and copied RFLO_ModMoveGridNconform to RFLO_ModMoveGridFrame
!
! Revision 1.8  2005/06/14 11:20:35  wasistho
! set itype to 1
!
! Revision 1.7  2005/06/13 21:47:35  wasistho
! changed patch%bcCoupled to patch%bcMotion
!
! Revision 1.6  2005/06/13 01:18:57  wasistho
! bug fixed in MgFrameEdge, ic to grid%ijkCorn(ic)
!
! Revision 1.5  2005/06/12 10:56:36  wasistho
! uncommented second/end TFI procedure
!
! Revision 1.4  2005/06/12 10:12:01  wasistho
! commented second/end TFI procedure
!
! Revision 1.3  2005/06/12 07:59:36  wasistho
! working version of Nconform1
!
! Revision 1.19  2005/06/12 06:21:29  wasistho
! modified dumax in MgFrameCorrectNeighbors
!
! Revision 1.18  2005/06/12 00:50:11  wasistho
! fixed bug in defining regNc
!
! Revision 1.17  2005/06/11 22:57:21  wasistho
! working version of RFLO_ModMoveGridNconform
!
! Revision 1.1  2005/06/11 20:59:26  wasistho
! import as non-conforming partition version of RFLO_ModMoveGridFrame
!
! Revision 1.13  2005/06/06 23:04:59  wasistho
! fixed bug put logical test between brackets
!
! Revision 1.12  2005/06/05 23:03:25  wasistho
! distinguish external boundary to be fully and partly external
!
! Revision 1.11  2005/06/04 22:28:06  wasistho
! more rigorous searching and made alghorithm 1 as default
!
! Revision 1.10  2005/06/04 10:51:43  wasistho
! move surface containing partly external-bc
!
! Revision 1.9  2005/06/04 04:25:44  wasistho
! fixed bug in determining interior vs external corners
!
! Revision 1.8  2005/06/04 01:01:33  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.7  2005/06/03 01:26:40  wasistho
! moved amplif to affect regCornOld gradually with iteration
!
! Revision 1.6  2005/06/02 22:59:34  wasistho
! added user controlled moveGridAmplif and moveGridPower
!
! Revision 1.5  2005/06/02 19:54:49  wasistho
! added control parameters amplif and pow
!
! Revision 1.4  2005/06/01 22:57:33  wasistho
! commented remeshing
!
! Revision 1.3  2005/06/01 08:02:39  wasistho
! increased mgFrame iteration from 5 to 10
!
! Revision 1.2  2005/06/01 07:14:41  wasistho
! debuged and made more robust
!
! Revision 1.1  2005/05/28 08:11:56  wasistho
! import RFLO_ModMoveGridFrame
!
!
!
! ******************************************************************************




















