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
! Purpose: provide initial solution for GENX.
!
! Description: the initialization is done separately for each region.
!              Flow conditions can vary between the regions.
!
! Input: gridLevel  = initial grid level
!        regions    = data for all regions
!        iReg       = region number
!        wins, winv = surface and volume window names
!
! Output: surf- and vol-grid, and initial cv.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_GenxInitSolution.F90,v 1.6 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GenxInitSolution( gridLevel,iReg,regions,wins,winv )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE PREP_ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
      RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, RFLO_ReadGridRegion, &
      RFLO_InitGenxInterfacePrep, InitializeFlowField
  USE ModParameters
  IMPLICIT NONE

  INCLUDE "comf90.h"
#include "Indexing.h"

! ... parameters
  INTEGER :: gridLevel, iReg
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkN, ng1, ng2

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch
  CHARACTER(CHRLEN)       :: wins, winv

  INTEGER :: iLev, lbound, iNOff, ijNOff
  INTEGER :: ibn, ien, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, dims(2), errorFlag

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'GenxInitSolution',&
  'PREP_GenxInitSolution.F90' )

! read grid and initialize solution -------------------------------------------

  global%winName = "Rocflo"

  wins = TRIM(global%winName)//'_surf'
  winv = TRIM(global%winName)//'_vol'

! Create volume and surface windows
  CALL COM_new_window( TRIM(wins) )
  CALL COM_new_window( TRIM(winv) )
  
! allocate memory for grid and read volume grid (already done by PREP_GetGrid)

  iLev = gridLevel

!  CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
!                                 jdnbeg,jdnend,kdnbeg,kdnend )
!  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
!  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
!  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
!  ALLOCATE( regions(iReg)%levels(iLev)%grid%xyz(3,ibn:ien),stat=errorFlag )
!  global%error = errorFlag
!  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

!  CALL RFLO_ReadGridRegion( iReg,regions )

! define surface grid

  DO iPatch=1,regions(iReg)%nPatches

    patch   => regions(iReg)%levels(iLev)%patches(iPatch)
!    lbound  =  patch%lbound

! - get dimensions, allocate memory

    ALLOCATE( patch%bcFlag(1),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

!    dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
!    dims(2) = ABS(patch%l2end-patch%l2beg) + 2
!    ALLOCATE( patch%surfCoord(3,dims(1),dims(2)),stat=errorFlag )
!    global%error = errorFlag
!    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

!    CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
!                                    ibeg,iend,jbeg,jend,kbeg,kend )

! - copy coordinates to temporary array

!    DO k=kbeg,kend
!      DO j=jbeg,jend
!        DO i=ibeg,iend
!          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
!          IF      (lbound==1 .OR. lbound==2) THEN
!            IF (lbound == 2) THEN
!              ng1 = j - jbeg + 1
!            ELSE
!              ng1 = jend - j + 1
!            ENDIF
!            ng2 = k - kbeg + 1
!          ELSE IF (lbound==3 .OR. lbound==4) THEN
!            ng1 = k - kbeg + 1
!            IF (lbound == 4) THEN
!              ng2 = i - ibeg + 1
!            ELSE
!              ng2 = iend - i + 1
!            ENDIF
!          ELSE IF (lbound==5 .OR. lbound==6) THEN
!            IF (lbound == 6) THEN
!              ng1 = i - ibeg + 1
!            ELSE
!              ng1 = iend - i + 1
!            ENDIF
!            ng2 = j - jbeg + 1
!          ENDIF
!          patch%surfCoord(1,ng1,ng2) = &
!                    regions(iReg)%levels(iLev)%grid%xyz(XCOORD,ijkN)
!          patch%surfCoord(2,ng1,ng2) = &
!                    regions(iReg)%levels(iLev)%grid%xyz(YCOORD,ijkN)
!          patch%surfCoord(3,ng1,ng2) = &
!                    regions(iReg)%levels(iLev)%grid%xyz(ZCOORD,ijkN)
!        ENDDO  ! i
!      ENDDO    ! j
!    ENDDO      ! k
  ENDDO        ! iPatch

! generating solution for Genx

  CALL InitializeFlowField( gridLevel,regions(iReg) )

! init. interface; restore grid and flow solution if restart ----------------

  WRITE(STDOUT,'(A)') SOLVER_NAME//' Preparing GenX interface ...'

  CALL RFLO_InitGenxInterfacePrep( iReg,regions(iReg),wins,winv )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GenxInitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_GenxInitSolution.F90,v $
! Revision 1.6  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/05/02 20:31:28  wasistho
! bug fixed, activated iLev
!
! Revision 1.3  2005/05/02 18:08:13  wasistho
! added cylindrical Taylor inflow profile capability
!
! Revision 1.2  2004/12/03 03:28:36  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.3  2004/11/13 22:59:28  wasistho
! invert orientation of genx-surface-variables
!
! Revision 1.2  2004/06/30 04:07:26  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.1  2004/06/30 00:06:05  wasistho
! initial import for GEN3
!
!
!******************************************************************************







