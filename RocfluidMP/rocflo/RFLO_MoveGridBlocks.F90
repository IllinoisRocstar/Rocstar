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
!          boundaries. This function distributes the deformations
!          in block-wise manner without any global smoothing.
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
! $Id: RFLO_MoveGridBlocks.F90,v 1.6 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridBlocks( regions )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY       : t_grid
  USE ModInterfaces, ONLY : RFLO_MoveGridSurfaces, RFLO_MoveGridInterfaces, &
        RFLO_ChangeInteriorGrid, RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
        RFLO_CalcCellCentroids, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
        RFLO_CheckMetrics, RFLO_CalcGridSpeeds, RFLO_GridRemesh
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
  INTEGER :: iReg

! ... local variables
  LOGICAL :: someMoved, someRemesh

  INTEGER :: iRemesh, jRemesh, nRemesh

  TYPE(t_grid), POINTER   :: grid, gridOld
  TYPE(t_global), POINTER :: global
#ifdef GENX
  DOUBLE PRECISION :: dAlpha
#endif

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MoveGridBlocks',&
  'RFLO_MoveGridBlocks.F90' )

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

! recompute grid metrics due to remeshing -------------------------------------

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

END SUBROUTINE RFLO_MoveGridBlocks

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MoveGridBlocks.F90,v $
! Revision 1.6  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/05 19:13:46  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.3  2005/05/28 05:52:38  wasistho
! refrain from grid remeshing
!
! Revision 1.2  2005/05/27 01:51:24  wasistho
! added remeshing
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.8  2004/09/02 02:35:20  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.7  2004/08/03 22:46:27  wasistho
! added RFLO_c2eAvgCoeffs
!
! Revision 1.6  2004/07/30 17:26:18  wasistho
! provide cell2face averaging coefficients
!
! Revision 1.5  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.15  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.14  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.13  2003/05/06 20:05:39  jblazek
! Corrected bug in grid motion (corner "averaging").
!
! Revision 1.12  2003/04/10 00:30:06  jblazek
! Corrected bug in grid movement algorithm.
!
! Revision 1.11  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.10  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.9  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.8  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/08/29 23:32:28  jblazek
! Changed name of index pointer for grid speeds.
!
! Revision 1.6  2002/08/27 20:44:31  jblazek
! Implemented calculation of grid speeds.
!
! Revision 1.5  2002/08/20 00:13:24  jblazek
! Removed print statements.
!
! Revision 1.4  2002/08/20 00:07:55  jblazek
! Corrected dimension of arcLen56.
!
! Revision 1.3  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







