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
! Purpose: set values in edge and corner cells of a region for smoke.
!
! Description: this is for the case if the cells are not located
!              within another region.
!
! Input: region = region data.
!
! Output: region%levels%peul%cv = cv variables in dummy cells.
!
! Notes: values in ALL edge and corner cells are averaged from neighboring
!        cells. Hence, this routine should be called first, before values
!        from adjacent regions are copied.
!
!******************************************************************************
!
! $Id: PEUL_SetCornerEdgeCells.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SetCornerEdgeCells( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijkD, ijkC1, ijkC2

  REAL(RFREAL),   POINTER :: cv(:,:)
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_SetCornerEdgeCells.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_SetCornerEdgeCells',&
  'PEUL_SetCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv => region%levels(iLev)%peul%cv

! edges 9, 10, 11, 12 ---------------------------------------------------------

  DO i=ipcbeg,ipcend
    DO j=jpcbeg-1,jdcbeg,-1         ! edges 9, 10
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i,j  ,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j+1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k+1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j+1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k-1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO

    DO j=jpcend+1,jdcend           ! edges 11, 12
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k+1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k-1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO
  ENDDO

! edges 2, 4, 6, 8 ------------------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg-1,idcbeg,-1        ! edges 2, 4
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k-1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k+1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO

    DO i=ipcend+1,idcend           ! edges 6, 8
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k-1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k+1,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO
  ENDDO

! edges 1, 3, 5, 7 ------------------------------------------------------------

  DO k=kdcbeg,kdcend
    DO i=ipcbeg-1,idcbeg,-1        ! edges 1, 3
      DO j=jpcbeg-1,jdcbeg,-1
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j+1,k,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO j=jpcend+1,jdcend
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j-1,k,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO

    DO i=ipcend+1,idcend           ! edges 5, 7
      DO j=jpcbeg-1,jdcbeg,-1
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j+1,k,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
      DO j=jpcend+1,jdcend
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j-1,k,iCOff,ijCOff)
        cv(:,ijkD) = 0.5_RFREAL*(cv(:,ijkC1) + cv(:,ijkC2))
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SetCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SetCornerEdgeCells.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:58  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/04/15 16:04:03  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.1  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
!******************************************************************************







