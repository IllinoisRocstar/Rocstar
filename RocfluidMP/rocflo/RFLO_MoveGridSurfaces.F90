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
!
! $Id: RFLO_MoveGridSurfaces.F90,v 1.5 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridSurfaces( regions,someMoved )

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

  CALL RegisterFunction( global,'RFLO_MoveGridSurfaces',&
  'RFLO_MoveGridSurfaces.F90' )

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

      CALL RFLO_EdgeDeformation( regions(iReg),grid%boundMoved,grid%edgeMoved, &
                                 grid%arcLen12,grid%arcLen34,grid%arcLen56, &
                                 gridOld%xyzOld,grid%xyz )

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

END SUBROUTINE RFLO_MoveGridSurfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MoveGridSurfaces.F90,v $
! Revision 1.5  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/14 04:38:49  wasistho
! added RFLO_EdgeDeformationStraight
!
! Revision 1.2  2006/03/05 19:02:53  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
!******************************************************************************







