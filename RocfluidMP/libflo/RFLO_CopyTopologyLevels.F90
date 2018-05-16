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
! Purpose: copy topology and BC data to all coarse grids.
!
! Description: none.
!
! Input: regions = dimensions and topology on finest grid.
!
! Output: regions = dimensions and topology on all coarse grids.
!
! Notes: conducted only for currently active regions on the own processor.
!
!******************************************************************************
!
! $Id: RFLO_CopyTopologyLevels.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CopyTopologyLevels( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_CopyBoundaryData
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iPatch

! ... local variables
  INTEGER :: nLevels, nPatches, ipc, jpc, kpc, bcType, errorFlag

  TYPE(t_patch), POINTER  :: patch, patchPrev
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_CopyTopologyLevels',&
  'RFLO_CopyTopologyLevels.F90' )

! loop over all regions

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      nLevels  = regions(iReg)%nGridLevels
      nPatches = regions(iReg)%nPatches

! --- loop over coarse grid levels

      DO iLev=2,nLevels
        ipc = regions(iReg)%levels(iLev-1)%grid%ipc     ! no. of physical cells
        jpc = regions(iReg)%levels(iLev-1)%grid%jpc
        kpc = regions(iReg)%levels(iLev-1)%grid%kpc
        regions(iReg)%levels(iLev)%grid%ipc = ipc/2
        regions(iReg)%levels(iLev)%grid%jpc = jpc/2
        regions(iReg)%levels(iLev)%grid%kpc = kpc/2

! ----- loop over patches

        ALLOCATE( regions(iReg)%levels(iLev)%patches(nPatches),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

        DO iPatch=1,nPatches
          patchPrev => regions(iReg)%levels(iLev-1)%patches(iPatch)
          patch     => regions(iReg)%levels(iLev  )%patches(iPatch)

! ------- current patch

          patch%bcType    = patchPrev%bcType
          patch%bcCoupled = patchPrev%bcCoupled
          patch%srcRegion = patchPrev%srcRegion
          patch%srcPatch  = patchPrev%srcPatch
          patch%align     = patchPrev%align
          patch%lbound    = patchPrev%lbound
          patch%l1beg     = (patchPrev%l1beg-1)/2 + 1
          patch%l1end     = patchPrev%l1end/2
          patch%l2beg     = (patchPrev%l2beg-1)/2 + 1
          patch%l2end     = patchPrev%l2end/2
          patch%srcLbound = patchPrev%srcLbound

          patch%mixt%bcSet     = patchPrev%mixt%bcSet
          patch%mixt%distrib   = patchPrev%mixt%distrib
          patch%mixt%nData     = patchPrev%mixt%nData
          patch%mixt%nSwitches = patchPrev%mixt%nSwitches

! ------- source patch (if there is an adjacent region)

          bcType = patch%bcType

          IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR.&
              (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR.&
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR.&
              (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR.&
              (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
            IF (patch%srcL1beg < patch%srcL1end) THEN
              patch%srcL1beg = (patchPrev%srcL1beg-1)/2 + 1
              patch%srcL1end = patchPrev%srcL1end/2
            ELSE
              patch%srcL1beg = patchPrev%srcL1beg/2
              patch%srcL1end = (patchPrev%srcL1end-1)/2 + 1
            ENDIF
            IF (patch%srcL2beg < patch%srcL2end) THEN
              patch%srcL2beg = (patchPrev%srcL2beg-1)/2 + 1
              patch%srcL2end = patchPrev%srcL2end/2
            ELSE
              patch%srcL2beg = patchPrev%srcL2beg/2
              patch%srcL2end = (patchPrev%srcL2end-1)/2 + 1
            ENDIF
          ELSE                                     ! no adjacent region
            patch%srcL1beg = patchPrev%srcL1beg
            patch%srcL1end = patchPrev%srcL1end
            patch%srcL2beg = patchPrev%srcL2beg
            patch%srcL2end = patchPrev%srcL2end
          ENDIF  ! bcType

! ------- distributions of boundary values

          IF ((bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR.&
              (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR.&
              (bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR.&
              (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR.&
              (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR.&
              (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE)) THEN
            CALL RFLO_CopyBoundaryData( global,patchPrev,patch )
          ENDIF  ! bcType

        ENDDO    ! iPatch

      ENDDO      ! iLev

    ENDIF        ! active region on my processor
  ENDDO          ! iReg

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_CopyTopologyLevels

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CopyTopologyLevels.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:38:00  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2002/10/25 18:36:19  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.6  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.5  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.4  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.3  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
! Revision 1.2  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







