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
! Purpose: update values in dummy cells for smoke at outflow boundary patch.
!
! Description: none.
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
! $Id: PEUL_BcondOutflow.F90,v 1.5 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_BcondOutflow( region,patch )

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
  INTEGER :: iLev, lbound, iCOff, ijCOff, ijkC, ijkC1, ijkD

  REAL(RFREAL), POINTER :: cv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_BcondOutflow.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_BcondOutflow',&
  'PEUL_BcondOutflow.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv => region%levels(iLev)%peul%cv

! loop over all cells of a patch

  DO idum=1,region%nDumCells

    IF (idum==1) THEN
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkD  = &
            IndIJK(i- idum   *idir,j- idum   *jdir,k- idum   *kdir,iCOff,ijCOff)
            ijkC  = &
            IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = &
            IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            cv(:,ijkD) = 2._RFREAL*cv(:,ijkC ) - cv(:,ijkC1) ! to be changed by
                                                             ! convective bc
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ELSE
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkD  = &
            IndIJK(i- idum   *idir,j- idum   *jdir,k- idum   *kdir,iCOff,ijCOff)
            ijkC  = &
            IndIJK(i-(idum-1)*idir,j-(idum-1)*jdir,k-(idum-1)*kdir,iCOff,ijCOff)
            ijkC1 = &
            IndIJK(i-(idum-2)*idir,j-(idum-2)*jdir,k-(idum-2)*kdir,iCOff,ijCOff)

            cv(:,ijkD) = 2._RFREAL*cv(:,ijkC ) - cv(:,ijkC1)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! idum
  ENDDO          ! idum

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_BcondOutflow

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_BcondOutflow.F90,v $
! Revision 1.5  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/16 20:36:50  wasistho
! reset to linear extrap. for now
!
! Revision 1.2  2006/03/15 19:50:30  wasistho
! improve peul outflow
!
! Revision 1.1  2004/12/01 21:09:17  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/09 15:12:04  jferry
! miscellaneous stylistic changes
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







