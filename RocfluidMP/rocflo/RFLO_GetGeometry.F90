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
! Purpose: read grid from file, define geometry for dummy cells.
!
! Description: none.
!
! Input: regions = grid dimensions
!        input from file.
!
! Output: regions%grid%xyz = grid coordinates.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_GetGeometry.F90,v 1.22 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetGeometry( regions,iread )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModInterfaces, ONLY : RFLO_ReadGrid, RFLO_GenerateCoarseGrids, &
                   RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
                   RFLO_ExchangeGeometry

  USE RFLO_ModGridMetrics,    ONLY : RFLO_ArcLengthPatch
  USE RFLO_ModGridRegionShape, ONLY: RFLO_GridFlatPatch, &
                                     RFLO_FindFunkyBlocks
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridControlMap3D, &
                                     RFLO_GridControlMap2D, &
                                     RFLO_GridControlGrad3D, &
                                     RFLO_GridControlGrad2D, &
                                     RFLO_GridControlFunc3D, &
                                     RFLO_GridControlFunc2D
  USE RFLO_ModForcesMoments,  ONLY : RFLO_FindPatchCoeffsGlo, &
                                     RFLO_WritePatchCoeffsInfo 
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)
  INTEGER :: iread

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  LOGICAL :: moveGrid
  REAL(RFREAL), POINTER :: xyzRef(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_GetGeometry',&
  'RFLO_GetGeometry.F90' )

! read/receive grid data for all regions if grid motion active

  moveGrid = .false.
  DO iReg=1,global%nRegions
    IF (regions(iReg)%mixtInput%moveGrid) moveGrid = .true.
  ENDDO

  IF (iread==1) THEN
#ifdef GENX
    IF (moveGrid) CALL RFLO_ReadGrid( regions )
#else
    CALL RFLO_ReadGrid( regions )
#endif
    CALL RFLO_FindFunkyBlocks( regions )
  ENDIF

! loop over all regions, generate coordinates of dummy nodes

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL RFLO_GenerateCoarseGrids( regions(iReg) )    ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )      ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )    ! extrapolate
    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

! exchange geometry between regions

  CALL RFLO_ExchangeGeometry( regions )

! store initial grid in separate arrays

  IF (iread==1 .AND. moveGrid) THEN
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor
        IF (regions(iReg)%mixtInput%moveGrid) THEN

          xyzRef => regions(iReg)%levels(1)%gridOld%xyzOld
          xyzRef =  regions(iReg)%levels(1)%grid%xyz
          
          regions(iReg)%levels(1)%grid%boundFlat(:)    = .TRUE.
          regions(iReg)%levels(1)%grid%edgeStraight(:) = .TRUE.
 
          DO iPatch=1,regions(iReg)%nPatches
            patch  => regions(iReg)%levels(1)%patches(iPatch)
            CALL RFLO_GridFlatPatch( regions(iReg),patch )
            CALL RFLO_ArcLengthPatch( regions(iReg),patch,xyzRef )
          ENDDO  ! iPatch
 
          IF (global%moveGridScheme==MOVEGRID_ELGLOBAL .OR. &
              global%moveGridScheme==MOVEGRID_ELFRAME) THEN

            CALL RFLO_GridControlMap3D( regions(iReg) )
            CALL RFLO_GridControlGrad3D( regions(iReg) )
            CALL RFLO_GridControlFunc3D( regions(iReg) )

            DO iPatch=1,regions(iReg)%nPatches
              patch  => regions(iReg)%levels(1)%patches(iPatch)
              IF (patch%bndFlat) THEN
                CALL RFLO_GridControlMap2D( regions(iReg),patch,iPatch )
                CALL RFLO_GridControlGrad2D( regions(iReg),patch,iPatch )
                CALL RFLO_GridControlFunc2D( regions(iReg),patch,iPatch )
              ENDIF ! flatPatch
            ENDDO   ! iPatch
          ENDIF     ! elliptic PDE

        ENDIF   ! moveGrid
      ENDIF     ! region on this processor and active
    ENDDO       ! iReg
  ENDIF         ! iread

! identify regions and patches contributing to global aero coeffs.
! based on original geometry and write to file

  IF (iread==1 .AND. global%aeroCoeffs==ACTIVE) THEN
    CALL RFLO_FindPatchCoeffsGlo( regions )
    CALL RFLO_WritePatchCoeffsInfo( regions )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetGeometry

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetGeometry.F90,v $
! Revision 1.22  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.21  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.20  2006/05/03 08:32:02  wasistho
! added if (patch%bndFlat) condition for 2D grid control
!
! Revision 1.19  2006/03/24 23:30:13  wasistho
! added FindPatchCoeffs and WritePatchCoeffs..
!
! Revision 1.18  2006/03/18 11:07:14  wasistho
! moved gridFlatPatch and findFunky.. to ModGridRegionShape
!
! Revision 1.17  2006/03/18 08:20:22  wasistho
! called ArcLengthPatch
!
! Revision 1.16  2006/03/14 04:34:04  wasistho
! initialized edgeStraight
!
! Revision 1.15  2006/03/12 22:07:02  wasistho
! changed RFLO_ElliptFlatPatch to RFLO_GridFlatPatch
!
! Revision 1.14  2006/03/12 10:29:54  wasistho
! initialized boundFlat
!
! Revision 1.13  2006/03/08 23:16:33  wasistho
! made gridOld%xyzOld available for all type gm
!
! Revision 1.12  2006/03/04 04:30:01  wasistho
! added RFLO_FindFunkyBlocks
!
! Revision 1.11  2006/03/02 01:27:32  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.10  2006/03/02 00:23:22  wasistho
! prepared elliptic pde grid motion
!
! Revision 1.9  2006/02/11 03:35:55  wasistho
! added calls controlGrad/Func 3D and 2D
!
! Revision 1.8  2006/02/08 07:52:15  wasistho
! added iPatch in controlMap2 argument
!
! Revision 1.7  2006/01/20 08:45:40  wasistho
! read .grda in Genx only for moving grid
!
! Revision 1.6  2005/12/07 08:47:14  wasistho
! added calls for surface mesh motion EPDE
!
! Revision 1.5  2005/12/05 10:48:07  wasistho
! added call RFLO_ElliptFlatPatch
!
! Revision 1.4  2005/12/03 09:34:11  wasistho
! compute control functions for movegrid EPDE
!
! Revision 1.3  2005/11/28 20:04:48  wasistho
! assigned gridOld%xyzOld
!
! Revision 1.2  2005/05/27 08:08:27  wasistho
! allow genx read initial grid
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.10  2004/07/18 21:54:28  jiao
! Updated not to call RFLO_ReadGrid in GENX.
!
! Revision 1.9  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.1  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







