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
! Purpose: correct values in the corner and edge cells at symmetry boundaries
!          for smoke.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch
!        bcType = type of boundary condition.
!
! Output: region%levels%peul%cv = cv variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_CorrectCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_CorrectCornerEdgeCells( region,patch,bcType )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
        RFLO_GetPatchIndices, RFLO_GetPatchDirection, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_patch),  INTENT(IN)    :: patch
  INTEGER,        INTENT(IN)    :: bcType

! ... loop variables
  INTEGER :: idum, i, j, k, iside

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  LOGICAL :: doDummy(4)

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ibegp, iendp, jbegp, jendp, kbegp, kendp, idir, jdir, kdir
  INTEGER :: ibeg(4), iend(4), jbeg(4), jend(4), kbeg(4), kend(4), nsides
  INTEGER :: iLev, lbound, iCOff, ijCOff
  INTEGER :: ijkC, ijkD, idm1

  REAL(RFREAL),   POINTER :: cv(:,:)
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_CorrectCornerEdgeCells.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_CorrectCornerEdgeCells',&
  'PEUL_CorrectCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibegp,iendp,jbegp,jendp,kbegp,kendp )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv => region%levels(iLev)%peul%cv

! set loop indices ------------------------------------------------------------

  doDummy(:) = .false.
  nsides     = 0

  IF      (lbound==1 .OR. lbound==2) THEN
    ibeg(:) = ibegp
    iend(:) = iendp
    IF (jbegp == jpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jdcbeg
      jend(nsides)    = jpcbeg - 1
      IF (kbegp == kpcbeg) THEN
        kbeg(nsides)  = kdcbeg
      ELSE
        kbeg(nsides)  = kbegp
      ENDIF
      IF (kendp == kpcend) THEN
        kend(nsides)  = kdcend
      ELSE
        kend(nsides)  = kendp
      ENDIF
    ENDIF
    IF (jendp == jpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jpcend + 1
      jend(nsides)    = jdcend
      IF (kbegp == kpcbeg) THEN
        kbeg(nsides)  = kdcbeg
      ELSE
        kbeg(nsides)  = kbegp
      ENDIF
      IF (kendp == kpcend) THEN
        kend(nsides)  = kdcend
      ELSE
        kend(nsides)  = kendp
      ENDIF
    ENDIF
    IF (kbegp == kpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jbegp
      jend(nsides)    = jendp
      kbeg(nsides)    = kdcbeg
      kend(nsides)    = kpcbeg - 1
    ENDIF
    IF (kendp == kpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      jbeg(nsides)    = jbegp
      jend(nsides)    = jendp
      kbeg(nsides)    = kpcend + 1
      kend(nsides)    = kdcend
    ENDIF

  ELSE IF (lbound==3 .OR. lbound==4) THEN
    jbeg(:) = jbegp
    jend(:) = jendp
    IF (kbegp == kpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kdcbeg
      kend(nsides)    = kpcbeg - 1
      IF (ibegp == ipcbeg) THEN
        ibeg(nsides)  = idcbeg
      ELSE
        ibeg(nsides)  = ibegp
      ENDIF
      IF (iendp == ipcend) THEN
        iend(nsides)  = idcend
      ELSE
        iend(nsides)  = iendp
      ENDIF
    ENDIF
    IF (kendp == kpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kpcend + 1
      kend(nsides)    = kdcend
      IF (ibegp == ipcbeg) THEN
        ibeg(nsides)  = idcbeg
      ELSE
        ibeg(nsides)  = ibegp
      ENDIF
      IF (iendp == ipcend) THEN
        iend(nsides)  = idcend
      ELSE
        iend(nsides)  = iendp
      ENDIF
    ENDIF
    IF (ibegp == ipcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kbegp
      kend(nsides)    = kendp
      ibeg(nsides)    = idcbeg
      iend(nsides)    = ipcbeg - 1
    ENDIF
    IF (iendp == ipcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      kbeg(nsides)    = kbegp
      kend(nsides)    = kendp
      ibeg(nsides)    = ipcend + 1
      iend(nsides)    = idcend
    ENDIF

  ELSE IF (lbound==5 .OR. lbound==6) THEN
    kbeg(:) = kbegp
    kend(:) = kendp
    IF (ibegp == ipcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = idcbeg
      iend(nsides)    = ipcbeg - 1
      IF (jbegp == jpcbeg) THEN
        jbeg(nsides)  = jdcbeg
      ELSE
        jbeg(nsides)  = jbegp
      ENDIF
      IF (jendp == jpcend) THEN
        jend(nsides)  = jdcend
      ELSE
        jend(nsides)  = jendp
      ENDIF
    ENDIF
    IF (iendp == ipcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ipcend + 1
      iend(nsides)    = idcend
      IF (jbegp == jpcbeg) THEN
        jbeg(nsides)  = jdcbeg
      ELSE
        jbeg(nsides)  = jbegp
      ENDIF
      IF (jendp == jpcend) THEN
        jend(nsides)  = jdcend
      ELSE
        jend(nsides)  = jendp
      ENDIF
    ENDIF
    IF (jbegp == jpcbeg) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ibegp
      iend(nsides)    = iendp
      jbeg(nsides)    = jdcbeg
      jend(nsides)    = jpcbeg - 1
    ENDIF
    IF (jendp == jpcend) THEN
      nsides          = nsides + 1
      doDummy(nsides) = .true.
      ibeg(nsides)    = ibegp
      iend(nsides)    = iendp
      jbeg(nsides)    = jpcend + 1
      jend(nsides)    = jdcend
    ENDIF
  ENDIF

! loop over all cells of a patch ----------------------------------------------

  DO iside=1,nsides
    IF (doDummy(iside)) THEN

      DO idum=1,region%nDumCells
        idm1 = idum - 1
        DO k=kbeg(iside),kend(iside)
          DO j=jbeg(iside),jend(iside)
            DO i=ibeg(iside),iend(iside)
              ijkD = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
              ijkC = IndIJK(i+idm1*idir,j+idm1*jdir,k+idm1*kdir,iCOff,ijCOff)

! ----------- symmetry boundary

              IF (bcType>=BC_SYMMETRY .AND. &
                  bcType<=BC_SYMMETRY+BC_RANGE) THEN
                cv(:,ijkD) = cv(:,ijkC)
              ENDIF

            ENDDO  ! i
          ENDDO    ! j
        ENDDO      ! k
      ENDDO        ! idum

    ENDIF  ! doDummy
  ENDDO    ! iside

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_CorrectCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_CorrectCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:31  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.1  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
!******************************************************************************







