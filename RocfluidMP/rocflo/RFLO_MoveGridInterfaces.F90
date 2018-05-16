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
!
! $Id: RFLO_MoveGridInterfaces.F90,v 1.4 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridInterfaces( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY       : t_grid
  USE ModInterfaces, ONLY : RFLO_ExchangeDnodeCopy, RFLO_EdgeDeformation, &
        RFLO_BoundaryDeformation, RFLO_ExchangeDnodeSend, &
        RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests
  USE ModError
  USE ModParameters
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

  CALL RegisterFunction( global,'RFLO_MoveGridInterfaces',&
  'RFLO_MoveGridInterfaces.F90' )

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
              CALL RFLO_EdgeDeformation( regions(iReg),grid%boundMoved, &
                                         grid%edgeMoved,grid%arcLen12, &
                                         grid%arcLen34,grid%arcLen56, &
                                         gridOld%xyzOld,grid%xyz )
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
              CALL RFLO_EdgeDeformation( regions(iReg),grid%boundMoved, &
                                         grid%edgeMoved,grid%arcLen12, &
                                         grid%arcLen34,grid%arcLen56, &
                                         gridOld%xyzOld,grid%xyz )
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

END SUBROUTINE RFLO_MoveGridInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MoveGridInterfaces.F90,v $
! Revision 1.4  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/05 19:03:04  wasistho
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







