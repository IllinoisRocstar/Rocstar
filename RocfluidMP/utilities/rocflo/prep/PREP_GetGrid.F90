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
! Purpose: obtain volume and surface grid.
!
! Description: none.
!
! Input: gridLevel  = initial grid level
!        regions    = data for all regions
!        iReg       = region number
!
! Output: vol-grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_GetGrid.F90,v 1.7 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GetGrid( gridLevel,iReg,regions )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE PREP_ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
         RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, RFLO_ReadGridRegion, &
         RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
         RFLO_GenerateCoarseGrids
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: gridLevel,iReg
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkN, ng1, ng2

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

  INTEGER :: iLev, lbound, iNOff, ijNOff
  INTEGER :: ibn, ien, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, dims(2), errorFlag

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'GetGrid',&
  'PREP_GetGrid.F90' )

! read grid -------------------------------------------------------------------
  
! allocate memory for grid

  regions(iReg)%currLevel = gridlevel
  iLev                    = regions(iReg)%currLevel

  CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
  ALLOCATE( regions(iReg)%levels(iLev)%grid%xyz(3,ibn:ien),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! read volume grid

  CALL RFLO_ReadGridRegion( iReg,regions )

  CALL RFLO_CopyGeometryDummy( regions(iReg) )      ! copy to dummy nodes
  CALL RFLO_ExtrapolateGeometry( regions(iReg) )    ! extrapolate

! define surface grid

  DO iPatch=1,regions(iReg)%nPatches

    patch   => regions(iReg)%levels(iLev)%patches(iPatch)
    lbound  =  patch%lbound

! - get dimensions, allocate memory

    dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
    dims(2) = ABS(patch%l2end-patch%l2beg) + 2
    ALLOCATE( patch%surfCoord(3,dims(1),dims(2)),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                    ibeg,iend,jbeg,jend,kbeg,kend )

! - initialize bcSet

    patch%mixt%bcSet = .FALSE.

! - copy coordinates to temporary array

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
          IF      (lbound==1 .OR. lbound==2) THEN
            IF (lbound == 2) THEN
              ng1 = j - jbeg + 1
            ELSE
              ng1 = jend - j + 1
            ENDIF
            ng2 = k - kbeg + 1
          ELSE IF (lbound==3 .OR. lbound==4) THEN
            ng1 = k - kbeg + 1
            IF (lbound == 4) THEN
              ng2 = i - ibeg + 1
            ELSE
              ng2 = iend - i + 1
            ENDIF
          ELSE IF (lbound==5 .OR. lbound==6) THEN
            IF (lbound == 6) THEN
              ng1 = i - ibeg + 1
            ELSE
              ng1 = iend - i + 1
            ENDIF
            ng2 = j - jbeg + 1
          ENDIF
          patch%surfCoord(1,ng1,ng2) = &
                    regions(iReg)%levels(iLev)%grid%xyz(XCOORD,ijkN)
          patch%surfCoord(2,ng1,ng2) = &
                    regions(iReg)%levels(iLev)%grid%xyz(YCOORD,ijkN)
          patch%surfCoord(3,ng1,ng2) = &
                    regions(iReg)%levels(iLev)%grid%xyz(ZCOORD,ijkN)
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GetGrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_GetGrid.F90,v $
! Revision 1.7  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:41:06  mparmar
! Renamed patch variables
!
! Revision 1.4  2005/09/30 00:13:16  wasistho
! added copy geometry to dummy
!
! Revision 1.3  2005/05/03 08:17:30  wasistho
! initialized bcSet to false
!
! Revision 1.2  2005/05/03 03:20:17  wasistho
! enabled modified cyl.Taylor inflow profile
!
! Revision 1.1  2005/05/02 18:09:34  wasistho
! added cylindrical Taylor inflow profile capability
!
!
!
!******************************************************************************







