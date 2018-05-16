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
! Purpose: set a variable to zero in all dummy cells.
!
! Description: this is intended to prevent a change of the solution inside
!              the dummy cells (these are set by the boundary conditions).
!
! Input: var = variable.
!
! Output: var = zeroed out in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ZeroDummyCells.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ZeroDummyCells( region,var )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  REAL(RFREAL), POINTER :: var(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijk

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ZeroDummyCells',&
  'RFLO_ZeroDummyCells.F90' )

! get dimensions

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

! side 1 and 2

  DO i=idcbeg,ipcbeg-1
    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO
  DO i=ipcend+1,idcend
    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO

! side 3 and 4

  DO j=jdcbeg,jpcbeg-1
    DO k=kpcbeg,kpcend
      DO i=idcbeg,idcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO
  DO j=jpcend+1,jdcend
    DO k=kpcbeg,kpcend
      DO i=idcbeg,idcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO

! side 5 and 6

  DO k=kdcbeg,kpcbeg-1
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO
  DO k=kpcend+1,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijk        = IndIJK(i,j,k,iCOff,ijCOff)
        var(:,ijk) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ZeroDummyCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ZeroDummyCells.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
!******************************************************************************







