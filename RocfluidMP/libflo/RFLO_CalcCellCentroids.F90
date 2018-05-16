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
! Purpose: calculate cell centroids (optionally).
!
! Description: none.
!
! Input: region%levels%grid = dimensions, coordinates (current region)
!
! Output: region%levels%grid%cofg = cell centroids
!
! Notes: centroids are defined only for the physical cells.
!
!******************************************************************************
!
! $Id: RFLO_CalcCellCentroids.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcCellCentroids( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, CentroidHexa
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, corner(8), ijkCell

  REAL(RFREAL)          :: xyzHexa(3,8)
  REAL(RFREAL), POINTER :: xyz(:,:), cofg(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcCellCentroids',&
       'RFLO_CalcCellCentroids.F90' )

! loop over all grid levels

  DO iLev=1,region%nGridLevels

    IF (ASSOCIATED(region%levels(iLev)%grid%cofg)) THEN      ! CG needed
      CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
      CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

      xyz  => region%levels(iLev)%grid%xyz
      cofg => region%levels(iLev)%grid%cofg

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ijkCell   = IndIJK(i,j,k,iCOff,ijCOff)
            corner(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
            corner(2) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
            corner(3) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
            corner(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
            corner(5) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
            corner(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
            corner(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
            corner(8) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)

            xyzHexa(1:3,1) = xyz(1:3,corner(1))
            xyzHexa(1:3,2) = xyz(1:3,corner(2))
            xyzHexa(1:3,3) = xyz(1:3,corner(3))
            xyzHexa(1:3,4) = xyz(1:3,corner(4))
            xyzHexa(1:3,5) = xyz(1:3,corner(5))
            xyzHexa(1:3,6) = xyz(1:3,corner(6))
            xyzHexa(1:3,7) = xyz(1:3,corner(7))
            xyzHexa(1:3,8) = xyz(1:3,corner(8))

            CALL CentroidHexa( xyzHexa, &
                               cofg(XCOORD,ijkCell), &
                               cofg(YCOORD,ijkCell), &
                               cofg(ZCOORD,ijkCell) )
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDIF        ! cofg allocated?

  ENDDO          ! iLev

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcCellCentroids

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcCellCentroids.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:15  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.3  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************







