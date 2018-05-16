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
! Purpose: search for all patches of a region to be used
!          for computing the thrust.
!
! Description: none.
!
! Input: region%levels%grid = dimensions, grid data, face vectors.
!
! Output: patch%thrustCalc = patch marked for thruct calculation (T/F).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FindThrustPatches.F90,v 1.4 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_FindThrustPatches( region,iReg )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            RFLO_GetPatchIndices, RFLO_GetPatchDirection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: plane, bcType, iNOff, ijNOff, ijkN0, ijkN1

  LOGICAL :: found

  REAL(RFREAL), POINTER :: xyz(:,:)
  REAL(RFREAL) :: tol=1.0e-6

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_FindThrustPatches',&
  'RFLO_FindThrustPatches.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz   => region%levels(iLev)%grid%xyz
  plane =  global%thrustPlane

! thrust plane within region? use tolerance since plane will be at boundary ---

!  IF (global%thrustCoord<MINVAL(xyz(plane,:)) .OR. &
!      global%thrustCoord>MAXVAL(xyz(plane,:))) THEN
  IF (global%thrustCoord-MINVAL(xyz(plane,:)) < -tol .OR. &
      global%thrustCoord-MAXVAL(xyz(plane,:)) >  tol) THEN
    DO iPatch=1,region%nPatches
      region%levels(iLev)%patches(iPatch)%thrustCalc = .false.
    ENDDO
    GOTO 999
  ENDIF

! loop over patches -----------------------------------------------------------

  DO iPatch=1,region%nPatches

    patch => region%levels(iLev)%patches(iPatch)
    bcType = patch%bcType

    patch%thrustCalc = .false.

    IF ((bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW+   BC_RANGE) .OR. &
        (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT+ BC_RANGE) .OR. &
        (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN

! --- get dimensions of the patch

      lbound = patch%lbound

      CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                      jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

! --- loop over all nodes of the patch

      found = .true.
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN0 = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
            ijkN1 = IndIJK(i+idir,j+jdir,k+kdir,iNOff,ijNOff)
!            IF (MIN(xyz(plane,ijkN0),xyz(plane,ijkN1))>global%thrustCoord .OR.&
!                MAX(xyz(plane,ijkN0),xyz(plane,ijkN1))<global%thrustCoord) &
            IF (global%thrustCoord-MIN(xyz(plane,ijkN0),xyz(plane,ijkN1)) < -tol .OR.&
                global%thrustCoord-MAX(xyz(plane,ijkN0),xyz(plane,ijkN1)) >  tol) &
              found = .false.
          ENDDO
        ENDDO
      ENDDO
      IF (found) patch%thrustCalc = .true.

    ENDIF   ! bcType
    IF (global%verbLevel>=VERBOSE_HIGH .AND. patch%thrustCalc) &
      WRITE(STDOUT,'(A,I5,A,I3)') '   - region ',iReg,', patch ',iPatch

  ENDDO     ! iPatch

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_FindThrustPatches

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FindThrustPatches.F90,v $
! Revision 1.4  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/15 09:46:38  rfiedler
! Use a tolerance when checking whether nodes are on thrust plane.
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************







