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
! Purpose: check that the number of cells inside the physical domain is
!          the same or larger than the number of dummy cells.
!
! Description: none.
!
! Input: regions%levels%grid = coordinates (physical cells) and dimensions
!
! Output: regions%levels%grid%xyz = coordinates (dummy cells)
!
! Notes: only region interfaces with continuous grid are treated here.
!
!******************************************************************************
!
! $Id: RFLO_CheckMinimumCells.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckMinimumCells( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : 
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER :: bcType, iRegSrc, nDumCells, ipc, jpc, kpc

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_CheckMinimumCells',&
  'RFLO_CheckMinimumCells.F90' )

! loop over all regions

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- loop over all levels

      DO iLev=1,regions(iReg)%nGridLevels

! ----- loop over patches

        DO iPatch=1,regions(iReg)%nPatches
          patch  => regions(iReg)%levels(iLev)%patches(iPatch)
          bcType =  patch%bcType

! ------- inter-region boundary / periodicity

          IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR.&
              (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR.&
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR.&
              (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR.&
              (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
            iRegSrc   = patch%srcRegion
            ipc       = regions(iReg)%levels(iLev)%grid%ipc
            jpc       = regions(iReg)%levels(iLev)%grid%jpc
            kpc       = regions(iReg)%levels(iLev)%grid%kpc
            nDumCells = regions(iRegSrc)%nDumCells
            IF (nDumCells > MIN(ipc,jpc,kpc)) THEN  ! better check all dims.
              WRITE(msg,1000) iReg,iLev,iPatch,bcType
              CALL ErrorStop( global,ERR_NUMBER_CELLS,__LINE__,msg )
            ENDIF
          ENDIF
        ENDDO  ! iPatch

      ENDDO    ! iLev

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', level ',I1,', patch ',I3,', BC type ',I3)

END SUBROUTINE RFLO_CheckMinimumCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckMinimumCells.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/20 00:47:10  jblazek
! Purpose of the routine was missing ...
!
! Revision 1.3  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
!******************************************************************************







