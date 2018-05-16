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
! Purpose: exchange grid coordinates between regions or between periodic
!          boundaries.
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
! $Id: RFLO_ExchangeGeometry.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometry( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ExchangeGeometryCopy, &
        RFLO_ExchangeGeometrySend, RFLO_ExchangeGeometryRecv, &
        RFLO_ExchangeGeometryLevels, RFLO_ClearSendRequests, &
        RFLO_ExchangeGeometryPrepare
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: bcType, iRegSrc, iPatchSrc

  TYPE(t_patch), POINTER  :: patch, patchSrc
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ExchangeGeometry',&
  'RFLO_ExchangeGeometry.F90' )

! exchange orientation of l1/2-directions -------------------------------------

  CALL RFLO_ExchangeGeometryPrepare( regions )

! exchange/send geometry between/to patches -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- loop over patches (finest grid level)

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType

! ----- conforming inter-region boundary, translational
!       and rotational periodicity

        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)

          IF (regions(iRegSrc)%procid == global%myProcid) THEN
            CALL RFLO_ExchangeGeometryCopy( regions(iReg),regions(iRegSrc), &
                                            patch,patchSrc )
            CALL RFLO_ExchangeGeometryLevels( regions(iReg),iPatch )
          ELSE
            CALL RFLO_ExchangeGeometrySend( regions(iReg),regions(iRegSrc), &
                                            patch )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! receive geometry ------------------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- loop over patches (finest grid level)

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType

! ----- conforming inter-region boundary, translational
!       or rotational periodicity

        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
            CALL RFLO_ExchangeGeometryRecv( regions(iReg),regions(iRegSrc), &
                                            patch,patchSrc )
            CALL RFLO_ExchangeGeometryLevels( regions(iReg),iPatch )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! clear send requests ---------------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL RFLO_ClearSendRequests( regions,iReg,.true. )
    ENDIF
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ExchangeGeometry

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometry.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.17  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.12  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.11  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.10  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.9  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.8  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.7  2002/04/01 19:36:08  jblazek
! Added routine to clear send requests.
!
! Revision 1.6  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.3  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2002/01/02 16:04:20  jblazek
! Added routines to generate geometry for dummy cells.
!
! Revision 1.1  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
!******************************************************************************







