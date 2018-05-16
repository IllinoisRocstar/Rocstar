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
! Purpose: zero out movements obtained by LaplaceGridSolve along a boundary.
!
! Description: none.
!
! Input: region = data of current region, grid movements
!        patch  = current patch.
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_LaplaceGridPatch.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkNB

  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_LaplaceGridPatch',&
  'RFLO_LaplaceGridPatch.F90' )

! get dimensions and pointers

  iLev = 1

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  xyzOld => region%levels(iLev)%grid%xyzOld

! new = old

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
        xyz(XCOORD,ijkNB) = xyzOld(XCOORD,ijkNB)
        xyz(YCOORD,ijkNB) = xyzOld(YCOORD,ijkNB)
        xyz(ZCOORD,ijkNB) = xyzOld(ZCOORD,ijkNB)
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_LaplaceGridPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_LaplaceGridPatch.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/08/25 21:51:23  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.1  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
!******************************************************************************







