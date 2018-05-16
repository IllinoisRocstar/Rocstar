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
! Purpose: update values in dummy cells for smoke at symmetry boundary patch.
!
! Description: (1) smoke densities in dummy cells are set equal to values
!                  in the corresponding interior nodes
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%peul = smoke variables in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_BcondSlipWall.F90,v 1.4 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BcondSlipWall( region,patch )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_patch),  INTENT(IN)    :: patch

! ... loop variables
  INTEGER :: idum, i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, lbound, iCOff, ijCOff, ijkC, ijkC1, ijkD, dumExtrapol

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_BcondSlipWall.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_BcondSlipWall',&
  'PEUL_BcondSlipWall.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  dumExtrapol = patch%mixt%switches(BCSWI_SLIPW_EXTRAP)

  cv => region%levels(iLev)%peul%cv

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkC  = IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
          ijkC1 = IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)

          IF (dumExtrapol == EXTRAPOL_CONST) THEN
            cv(:,ijkD) = cv(:,ijkC)
          ELSE
            cv(:,ijkD) = 2._RFREAL*cv(:,ijkC ) - cv(:,ijkC1)
          END IF ! dumExtrapol

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BcondSlipWall

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BcondSlipWall.F90,v $
! Revision 1.4  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:19  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/12/01 21:09:18  haselbac
! Initial revision after changing case
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1  2003/04/09 15:09:56  jferry
! added slip wall boundary conditions
!
!******************************************************************************







