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
! Purpose: Suite of randvoorwarden.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModStatsBoundaryConditions.F90,v 1.4 2008/12/06 08:44:17 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModStatsBoundaryConditions

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModBndPatch, ONLY  : t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_StatBcondCopy, &
            RFLO_StatBcondNegate, &
            RFLO_StatBcondExtrap, &
            RFLO_StatBoundaryConditionsSet
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModStatsBoundaryConditions.F90,v $ $Revision: 1.4 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: update dummy cells values by copying from adjacent interior cells.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch
!        istbeg,istend = beginning and end statistics index being treated
!        tav = statistics variables in current region
!
! Output: tav = statistics variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_StatBcondCopy( region,patch,istbeg,istend,tav )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: istbeg,istend
  REAL(RFREAL), POINTER :: tav(:,:)

! ... loop variables
  INTEGER :: idum, i, j, k, ist

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, iCOff, ijCOff, lbound, ndumi, ndumj, ndumk
  INTEGER :: ijkC, ijkD

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_StatBcondCopy',&
  'RFLO_ModStatsBoundaryConditions.F90' )

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ndumi = region%nDumCells
  ndumj = region%nDumCells
  ndumk = region%nDumCells
  IF (lbound==1 .OR. lbound==2) ndumi = 0
  IF (lbound==3 .OR. lbound==4) ndumj = 0
  IF (lbound==5 .OR. lbound==6) ndumk = 0

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg-ndumk,kend+ndumk
      DO j=jbeg-ndumj,jend+ndumj
        DO i=ibeg-ndumi,iend+ndumi
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkC = IndIJK(i,j,k,iCOff,ijCOff)

          DO ist=istbeg,istend
            tav(ist,ijkD) = tav(ist,ijkC) 
          ENDDO

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_StatBcondCopy

!******************************************************************************
!
! Purpose: update dummy cells values by negating adjacent interior cells.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch
!        istbeg,istend = beginning and end statistics index being treated
!        tav = statistics variables in current region
!
! Output: tav = statistics variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_StatBcondNegate( region,patch,istbeg,istend,tav )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: istbeg,istend
  REAL(RFREAL), POINTER :: tav(:,:)

! ... loop variables
  INTEGER :: idum, i, j, k, ist

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, iCOff, ijCOff, lbound, ndumi, ndumj, ndumk
  INTEGER :: ijkC, ijkD

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_StatBcondNegate',&
  'RFLO_ModStatsBoundaryConditions.F90' )

! get dimensions and pointers

  iLev = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ndumi = region%nDumCells
  ndumj = region%nDumCells
  ndumk = region%nDumCells
  IF (lbound==1 .OR. lbound==2) ndumi = 0
  IF (lbound==3 .OR. lbound==4) ndumj = 0
  IF (lbound==5 .OR. lbound==6) ndumk = 0

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg-ndumk,kend+ndumk
      DO j=jbeg-ndumj,jend+ndumj
        DO i=ibeg-ndumi,iend+ndumi
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkC = IndIJK(i,j,k,iCOff,ijCOff)

          DO ist=istbeg,istend
            tav(ist,ijkD) = -tav(ist,ijkC) 
          ENDDO

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_StatBcondNegate

!******************************************************************************
!
! Purpose: update dummy cells values by extrapolation.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch
!        istbeg,istend = beginning and end statistics index being treated
!        tav = statistics variables in current region
!
! Output: tav = statistics variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_StatBcondExtrap( region,patch,istbeg,istend,tav )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: istbeg,istend
  REAL(RFREAL), POINTER :: tav(:,:)

! ... loop variables
  INTEGER :: idum, i, j, k, ist

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, iCOff, ijCOff, lbound, ndumi, ndumj, ndumk
  INTEGER :: ijkC, ijkC1, ijkD

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_StatBcondExtrap',&
  'RFLO_ModStatsBoundaryConditions.F90' )

! get dimensions and pointers

  iLev = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ndumi = region%nDumCells
  ndumj = region%nDumCells
  ndumk = region%nDumCells
  IF (lbound==1 .OR. lbound==2) ndumi = 0
  IF (lbound==3 .OR. lbound==4) ndumj = 0
  IF (lbound==5 .OR. lbound==6) ndumk = 0

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg-ndumk,kend+ndumk
      DO j=jbeg-ndumj,jend+ndumj
        DO i=ibeg-ndumi,iend+ndumi
          ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkC = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
          ijkC1= IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

          DO ist=istbeg,istend
            tav(ist,ijkD) = 2._RFREAL*tav(ist,ijkC)-tav(ist,ijkC1) 
          ENDDO

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_StatBcondExtrap

!******************************************************************************
!
! Purpose: set boundary conditions
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%mixt%tav, turb%tav, plag%tav, peul%tav = 
!         updated statistics of all physical modules if applicable in dummy 
!         cells of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_StatBoundaryConditionsSet( regions )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: iLev, nPatches, bcType, istbeg, istend, lbound
  REAL(RFREAL), POINTER :: tav(:,:)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_StatBoundaryConditionsSet',&
  'RFLO_ModStatsBoundaryConditions.F90' )

! start -----------------------------------------------------------------------

  IF (global%statBc==1) THEN

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor

! - get dimensions ------------------------------------------------------------

    iLev     = regions(iReg)%currLevel
    nPatches = regions(iReg)%nPatches

! - loop over patches ---------------------------------------------------------

    DO iPatch=1,nPatches

      patch => regions(iReg)%levels(iLev)%patches(iPatch)

      bcType = patch%bcType
      lbound = patch%lbound

! --- copy dummy

      IF ((bcType>=BC_INFLOW .AND. bcType<=BC_INFLOW+BC_RANGE) .OR. &
          (bcType>=BC_OUTFLOW .AND. bcType<=BC_OUTFLOW+BC_RANGE) .OR. &
          (bcType>=BC_FARFIELD .AND. bcType<=BC_FARFIELD+BC_RANGE) .OR. &
          (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE)) THEN

        IF (global%mixtNStat > 0) THEN
          tav => regions(iReg)%levels(iLev)%mixt%tav
          istbeg = 1
          istend = global%mixtNStat
          CALL RFLO_StatBcondCopy( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#ifdef TURB
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
            (global%turbNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%turb%tav
          istbeg = 1
          istend = global%turbNStat
          CALL RFLO_StatBcondCopy( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PLAG
        IF ((global%plagUsed .EQV. .TRUE.) .AND. &
            (global%plagNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%plag%tav
          istbeg = 1
          istend = global%plagNStat
          CALL RFLO_StatBcondCopy( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PEUL
!        IF ((global%peulUsed .EQV. .TRUE.) .AND. &
!            (global%peulNStat > 0)) THEN
!          tav => regions(iReg)%levels(iLev)%peul%tav
!          istbeg = 1
!          istend = global%peulNStat
!          CALL RFLO_StatBcondCopy( regions(iReg),patch,istbeg,istend,tav )
!        ENDIF
#endif

! --- slip wall

      ELSE IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) THEN

        IF (global%mixtNStat > 0) THEN
          tav => regions(iReg)%levels(iLev)%mixt%tav
          istbeg = 1
          istend = global%mixtNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          IF (lbound==1 .OR. lbound==2) THEN
            istbeg = 2   ! u
            istend = 2
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 7   ! uu
            istend = 7
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 10  ! uv
            istend = 10
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
        ENDIF
#ifdef TURB
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
            (global%turbNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%turb%tav
          istbeg = 1
          istend = global%turbNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PLAG
        IF ((global%plagUsed .EQV. .TRUE.) .AND. &
            (global%plagNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%plag%tav
          istbeg = 1
          istend = global%plagNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          IF (lbound==1 .OR. lbound==2) THEN
            istbeg = 4   ! u
            istend = 4
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 10  ! uu
            istend = 10
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
        ENDIF
#endif
#ifdef PEUL
!        IF ((global%peulUsed .EQV. .TRUE.) .AND. &
!            (global%peulNStat > 0)) THEN
!          tav => regions(iReg)%levels(iLev)%peul%tav
!          istbeg = 1
!          istend = global%peulNStat
!          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
!        ENDIF
#endif

! --- noslip wall

      ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

        IF (global%mixtNStat > 0) THEN
          tav => regions(iReg)%levels(iLev)%mixt%tav
          istbeg = 1
          istend = global%mixtNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          istbeg = 2   ! u,v,w
          istend = 4   
          CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          istbeg = 7   ! uu,vv,ww,uv
          istend = 10
          CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#ifdef TURB
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
            (global%turbNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%turb%tav
          istbeg = 1
          istend = global%turbNStat
          CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PLAG
        IF ((global%plagUsed .EQV. .TRUE.) .AND. &
            (global%plagNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%plag%tav
          istbeg = 1
          istend = global%plagNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          istbeg = 4   ! u,v,w
          istend = 6
          CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          istbeg = 10  ! uu
          istend = 10
          CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PEUL
!        IF ((global%peulUsed .EQV. .TRUE.) .AND. &
!            (global%peulNStat > 0)) THEN
!          tav => regions(iReg)%levels(iLev)%peul%tav
!          istbeg = 1
!          istend = global%peulNStat
!          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
!        ENDIF
#endif

! --- injection

      ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

        IF (global%mixtNStat > 0) THEN
          tav => regions(iReg)%levels(iLev)%mixt%tav
          istbeg = 1
          istend = global%mixtNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          IF (lbound==1 .OR. lbound==2) THEN
            istbeg = 3   ! v,w
            istend = 4
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 8   ! vv,ww,uv
            istend = 10
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
          IF (lbound/=1 .AND. lbound/=2) THEN
            istbeg = 2   ! u
            istend = 2
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 7   ! uu
            istend = 7
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 10  ! uv
            istend = 10
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
        ENDIF
#ifdef TURB
        IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
            (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) .AND. &
            (global%turbNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%turb%tav
          istbeg = 1
          istend = global%turbNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
        ENDIF
#endif
#ifdef PLAG
        IF ((global%plagUsed .EQV. .TRUE.) .AND. &
            (global%plagNStat > 0)) THEN
          tav => regions(iReg)%levels(iLev)%plag%tav
          istbeg = 1
          istend = global%plagNStat
          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
          IF (lbound==1 .OR. lbound==2) THEN
            istbeg = 5   ! v,w
            istend = 6
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
          IF (lbound/=1 .AND. lbound/=2) THEN
            istbeg = 4   ! u
            istend = 4
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
            istbeg = 10  ! uv
            istend = 10
            CALL RFLO_StatBcondNegate( regions(iReg),patch,istbeg,istend,tav )
          ENDIF
        ENDIF
#endif
#ifdef PEUL
!        IF ((global%peulUsed .EQV. .TRUE.) .AND. &
!            (global%peulNStat > 0)) THEN
!          tav => regions(iReg)%levels(iLev)%peul%tav
!          istbeg = 1
!          istend = global%peulNStat
!          CALL RFLO_StatBcondExtrap( regions(iReg),patch,istbeg,istend,tav )
!        ENDIF
#endif
      ENDIF  ! bcType

    ENDDO  ! iPatch
    ENDIF  ! myProcid
  ENDDO    ! iReg

  ENDIF    ! statBc

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_StatBoundaryConditionsSet

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLO_ModStatsBoundaryConditions

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModStatsBoundaryConditions.F90,v $
! Revision 1.4  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/06/16 22:43:11  wasistho
! get rid of bc_unknown option
!
! Revision 1.1  2005/05/21 01:41:54  wasistho
! added RFLO_ModStatsBc
!
!
! ******************************************************************************










