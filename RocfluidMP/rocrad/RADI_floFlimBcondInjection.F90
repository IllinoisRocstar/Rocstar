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
! Purpose: update radiant Energy in dummy cells at boundary patch assuming
!          injection boundary.
!
! Description: none.
!
! Input: region = region dimensions, user input
!        patch  = current patch.
!
! Output: region%levels%radi%cv = radiant Energy in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_floFlimBcondInjection.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimBcondInjection( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: idum, i, j, k, l

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iLev, iCOff, ijCOff, ijkC, ijkC1, ijkD, nCv

  REAL(RFREAL), POINTER :: cv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FloFlimBcondInjection',&
  'RADI_floFlimBcondInjection.F90' )

! get dimensions and pointers

  iLev   = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  nCv =  region%radiInput%nCv
  cv  => region%levels(iLev)%radi%cv

! loop over all cells of a patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkD  = IndIJK(i-idum*idir,j-idum*jdir,k-idum*kdir,iCOff,ijCOff)
          ijkC  = IndIJK(i, j, k,iCOff,ijCOff)
          ijkC1 = IndIJK(i+(idum-1)*idir,j+(idum-1)*jdir,k+(idum-1)*kdir,iCOff,ijCOff)
          IF (idum == 1) THEN
            DO l=1,nCv
              cv(l,ijkD) = cv(l,ijkC)
            ENDDO
          ELSE
            DO l=1,nCv
              cv(l,ijkD) = 2._RFREAL*cv(l,ijkC) - cv(l,ijkC1)
            ENDDO
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FloFlimBcondInjection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimBcondInjection.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







