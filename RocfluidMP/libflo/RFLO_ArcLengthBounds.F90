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
! Purpose: calculate approximate arclengths for every grid line
!          between two oposite region boundaries (on the finest
!          grid only).
!
! Description: none.
!
! Input: region = grid dimensions
!        xyz    = node coordinates.
!
! Output: arcLen12 = arclength between i=const. boundaries for each j, k
!         arcLen34 = arclength between j=const. boundaries for each k, i
!         arcLen56 = arclength between k=const. boundaries for each i, j.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ArcLengthBounds.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ArcLengthBounds( region,xyz,arcLen12,arcLen34,arcLen56 )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff, ijkN, ijkN1
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ArcLengthBounds',&
       'RFLO_ArcLengthBounds.F90' )

! get pointers and dimensions

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! zero out variables

  arcLen12(:,:) = 0._RFREAL
  arcLen34(:,:) = 0._RFREAL
  arcLen56(:,:) = 0._RFREAL

! boundaries 1 and 2 (i=const.)

  DO i=ipnbeg+1,ipnend
    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)
        ijkN1 = IndIJK(i-1,j,k,iNOff,ijNOff)
        arcLen12(j,k) = arcLen12(j,k) + &
                        SQRT((xyz(XCOORD,ijkN)-xyz(XCOORD,ijkN1))**2 + &
                             (xyz(YCOORD,ijkN)-xyz(YCOORD,ijkN1))**2 + &
                             (xyz(ZCOORD,ijkN)-xyz(ZCOORD,ijkN1))**2)
      ENDDO
    ENDDO
  ENDDO

! boundaries 3 and 4 (j=const.)

  DO j=jpnbeg+1,jpnend
    DO i=ipnbeg,ipnend
      DO k=kpnbeg,kpnend
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(i,j-1,k,iNOff,ijNOff)
        arcLen34(k,i) = arcLen34(k,i) + &
                        SQRT((xyz(XCOORD,ijkN)-xyz(XCOORD,ijkN1))**2 + &
                             (xyz(YCOORD,ijkN)-xyz(YCOORD,ijkN1))**2 + &
                             (xyz(ZCOORD,ijkN)-xyz(ZCOORD,ijkN1))**2)
      ENDDO
    ENDDO
  ENDDO

! boundarie 5 and 6 (k=const.)

  DO k=kpnbeg+1,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,j,k-1,iNOff,ijNOff)
        arcLen56(i,j) = arcLen56(i,j) + &
                        SQRT((xyz(XCOORD,ijkN)-xyz(XCOORD,ijkN1))**2 + &
                             (xyz(YCOORD,ijkN)-xyz(YCOORD,ijkN1))**2 + &
                             (xyz(ZCOORD,ijkN)-xyz(ZCOORD,ijkN1))**2)
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ArcLengthBounds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ArcLengthBounds.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:15  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2003/03/14 22:05:10  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:48:04  jblazek
! Implemented grid deformation capability.
!
!******************************************************************************







