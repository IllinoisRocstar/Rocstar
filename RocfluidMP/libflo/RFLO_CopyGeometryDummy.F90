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
! Purpose: copy grid coordinates from boundary to all dummy nodes.
!
! Description: none.
!
! Input: region%levels%grid = coordinates (physical cells) and dimensions
!
! Output: region%levels%grid%xyz = coordinates (dummy cells)
!
! Notes: boundaries between regions or periodic boundaries are treated
!        later. Thus, this routine does not distinguish between boundary
!        types.
!
!******************************************************************************
!
! $Id: RFLO_CopyGeometryDummy.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CopyGeometryDummy( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetDimensPhysNodes, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, ijk, ijkD

  REAL(RFREAL), POINTER :: xyz(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CopyGeometryDummy',&
  'RFLO_CopyGeometryDummy.F90' )

! loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    xyz => region%levels(iLev)%grid%xyz

! - side 1 and 2

    DO i=idnbeg,ipnbeg-1
      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          ijk  = IndIJK(ipnbeg,j,k,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i     ,j,k,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO
    DO i=ipnend+1,idnend
      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          ijk  = IndIJK(ipnend,j,k,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i     ,j,k,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO

! - side 3 and 4

    DO j=jdnbeg,jpnbeg-1
      DO k=kpnbeg,kpnend
        DO i=idnbeg,idnend
          ijk  = IndIJK(i,jpnbeg,k,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i,j     ,k,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO
    DO j=jpnend+1,jdnend
      DO k=kpnbeg,kpnend
        DO i=idnbeg,idnend
          ijk  = IndIJK(i,jpnend,k,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i,j     ,k,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO

! - side 5 and 6

    DO k=kdnbeg,kpnbeg-1
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend
          ijk  = IndIJK(i,j,kpnbeg,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i,j,k     ,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO
    DO k=kpnend+1,kdnend
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend
          ijk  = IndIJK(i,j,kpnend,iNOff,ijNOff)   ! boundary node
          ijkD = IndIJK(i,j,k     ,iNOff,ijNOff)   ! dummy node
          xyz(XCOORD,ijkD) = xyz(XCOORD,ijk)
          xyz(YCOORD,ijkD) = xyz(YCOORD,ijk)
          xyz(ZCOORD,ijkD) = xyz(ZCOORD,ijk)
        ENDDO
      ENDDO
    ENDDO

  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CopyGeometryDummy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CopyGeometryDummy.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
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
! Revision 1.3  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.2  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.1  2002/01/02 16:04:20  jblazek
! Added routines to generate geometry for dummy cells.
!
!******************************************************************************







