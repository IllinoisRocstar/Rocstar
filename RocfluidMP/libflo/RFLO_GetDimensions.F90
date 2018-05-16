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
! Purpose: calculate dimensions of the computational domain
!          (start & end index).
!
! Description: file contains the following subroutines:
!
!  - GetDimensDummy        = including dummy cells
!  - GetDimensPhys         = cells within physical domain only
!  - GetDimensDummyNodes   = including dummy nodes
!  - GetDimensPhysNodes    = grid nodes within physical domain only
!  - GetCellOffset         = total # of cells in i-direction and for combination
!                            of i- and j-direction (incl. dummy cells)
!  - GetNodeOffset         = total # of nodes in i-direction and for combination
!                            of i- and j-direction (incl. dummy nodes)
!  - GetEdgeCellsIndices   = dimensions of a "block" of edge cells
!  - GetCornerCellsIndices = dimensions of a "block" of corner cells
!
! Input: region  = current region
!        iLev    = current grid level for region
!        iedge   = edge of a hexahedron (number between 1 and 12)
!        icorner = corner of a hexahedron (number between 1 and 8)
!
! Output: idcbeg, ..., kdcend = dimensions including dummy cells
!         ipcbeg, ..., kpcend = dimensions including physical cells only
!         idnbeg, ..., kdnend = dimensions including dummy nodes
!         ipnbeg, ..., kpnend = dimensions including physical nodes only
!         iebeg, ..., keend   = start & end indices of edge cells
!         icbeg, ..., kcend   = start & end indices of corner cells
!         iCellOffset         = total number of cells in i-direction
!         ijCellOffset        = total number of cells in i*j-direction
!         iNodeOffset         = total number of nodes in i-direction
!         ijNodeOffset        = total number of nodes in i*j-direction
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_GetDimensions.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                                kdcbeg,kdcend )

  USE ModDataStruct, ONLY : t_region
  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  TYPE(t_region) :: region

! ... local variables
  CHARACTER(CHRLEN) :: msg

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetDimensDummy',&
  'RFLO_GetDimensions.F90' )

  IF (iLev<1 .OR. iLev>region%nGridLevels) THEN
    WRITE(msg,1000) iLev,region%nGridLevels
    CALL ErrorStop( region%global,ERR_GRID_LEVEL,&
    __LINE__,msg )
  ENDIF

  idcbeg = 1                            - region%nDumCells
  idcend = region%levels(iLev)%grid%ipc + region%nDumCells
  jdcbeg = 1                            - region%nDumCells
  jdcend = region%levels(iLev)%grid%jpc + region%nDumCells
  kdcbeg = 1                            - region%nDumCells
  kdcend = region%levels(iLev)%grid%kpc + region%nDumCells

  CALL DeregisterFunction( region%global )

1000 FORMAT('Grid level= ',I2,' but should be 1-',I1,'.')

END SUBROUTINE RFLO_GetDimensDummy

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                               kpcbeg,kpcend )

  USE ModDataStruct, ONLY : t_region
  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  TYPE(t_region) :: region

! ... local variables
  CHARACTER(CHRLEN) :: msg

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetDimensPhys',&
  'RFLO_GetDimensions.F90' )

  IF (iLev<1 .OR. iLev>region%nGridLevels) THEN
    WRITE(msg,1000) iLev,region%nGridLevels
    CALL ErrorStop( region%global,ERR_GRID_LEVEL,&
    __LINE__,msg )
  ENDIF

  ipcbeg = 1
  ipcend = region%levels(iLev)%grid%ipc
  jpcbeg = 1
  jpcend = region%levels(iLev)%grid%jpc
  kpcbeg = 1
  kpcend = region%levels(iLev)%grid%kpc

  CALL DeregisterFunction( region%global )

1000 FORMAT('Grid level= ',I2,' but should be 1-',I1,'.')

END SUBROUTINE RFLO_GetDimensPhys

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend, &
                                     kdnbeg,kdnend )

  USE ModDataStruct, ONLY : t_region
  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  TYPE(t_region) :: region

! ... local variables
  CHARACTER(CHRLEN) :: msg

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetDimensDummyNodes',&
  'RFLO_GetDimensions.F90' )

  IF (iLev<1 .OR. iLev>region%nGridLevels) THEN
    WRITE(msg,1000) iLev,region%nGridLevels
    CALL ErrorStop( region%global,ERR_GRID_LEVEL,&
    __LINE__,msg )
  ENDIF

  idnbeg = 1 - region%nDumCells
  idnend = 1 + region%nDumCells + region%levels(iLev)%grid%ipc
  jdnbeg = 1 - region%nDumCells
  jdnend = 1 + region%nDumCells + region%levels(iLev)%grid%jpc
  kdnbeg = 1 - region%nDumCells
  kdnend = 1 + region%nDumCells + region%levels(iLev)%grid%kpc

  CALL DeregisterFunction( region%global )

1000 FORMAT('Grid level= ',I2,' but should be 1-',I1,'.')

END SUBROUTINE RFLO_GetDimensDummyNodes

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                    kpnbeg,kpnend )

  USE ModDataStruct, ONLY : t_region
  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  TYPE(t_region) :: region

! ... local variables
  CHARACTER(CHRLEN) :: msg

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetDimensPhysNodes',&
  'RFLO_GetDimensions.F90' )

  IF (iLev<1 .OR. iLev>region%nGridLevels) THEN
    WRITE(msg,1000) iLev,region%nGridLevels
    CALL ErrorStop( region%global,ERR_GRID_LEVEL,&
    __LINE__,msg )
  ENDIF

  ipnbeg = 1
  ipnend = region%levels(iLev)%grid%ipc + 1
  jpnbeg = 1
  jpnend = region%levels(iLev)%grid%jpc + 1
  kpnbeg = 1
  kpnend = region%levels(iLev)%grid%kpc + 1

  CALL DeregisterFunction( region%global )

1000 FORMAT('Grid level= ',I2,' but should be 1-',I1,'.')

END SUBROUTINE RFLO_GetDimensPhysNodes

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetCellOffset( region,iLev,iCellOffset,ijCellOffset )

  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, iCellOffset, ijCellOffset
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetCellOffset',&
  'RFLO_GetDimensions.F90' )

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )

  iCellOffset  = idcend - idcbeg + 1
  ijCellOffset = iCellOffset*(jdcend-jdcbeg+1)

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetCellOffset

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )

  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, iNodeOffset, ijNodeOffset
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetNodeOffset',&
  'RFLO_GetDimensions.F90' )

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )

  iNodeOffset  = idnend - idnbeg + 1
  ijNodeOffset = iNodeOffset*(jdnend-jdnbeg+1)

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetNodeOffset

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetEdgeCellsIndices( region,iLev,iedge, &
                                     iebeg,ieend,jebeg,jeend,kebeg,keend )

  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, iedge
  INTEGER        :: iebeg, ieend, jebeg, jeend, kebeg, keend
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetEdgeCellsIndices',&
  'RFLO_GetDimensions.F90' )

  IF (iedge<0 .OR. iedge>12) &
    CALL ErrorStop( region%global,ERR_VOLUME_EDGES,&
    __LINE__ )

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )

  SELECT CASE (iedge)
    CASE (1)
      iebeg = idcbeg
      ieend = ipcbeg - 1
      jebeg = jdcbeg
      jeend = jpcbeg - 1
      kebeg = kpcbeg
      keend = kpcend
    CASE (2)
      iebeg = idcbeg
      ieend = ipcbeg - 1
      jebeg = jpcbeg
      jeend = jpcend
      kebeg = kpcend + 1
      keend = kdcend
    CASE (3)
      iebeg = idcbeg
      ieend = ipcbeg - 1
      jebeg = jpcend + 1
      jeend = jdcend
      kebeg = kpcbeg
      keend = kpcend
    CASE (4)
      iebeg = idcbeg
      ieend = ipcbeg - 1
      jebeg = jpcbeg
      jeend = jpcend
      kebeg = kdcbeg
      keend = kpcbeg - 1
    CASE (5)
      iebeg = ipcend + 1
      ieend = idcend
      jebeg = jdcbeg
      jeend = jpcbeg - 1
      kebeg = kpcbeg
      keend = kpcend
    CASE (6)
      iebeg = ipcend + 1
      ieend = idcend
      jebeg = jpcbeg
      jeend = jpcend
      kebeg = kpcend + 1
      keend = kdcend
    CASE (7)
      iebeg = ipcend + 1
      ieend = idcend
      jebeg = jpcend + 1
      jeend = jdcend
      kebeg = kpcbeg
      keend = kpcend
    CASE (8)
      iebeg = ipcend + 1
      ieend = idcend
      jebeg = jpcbeg
      jeend = jpcend
      kebeg = kdcbeg
      keend = kpcbeg - 1
    CASE (9)
      iebeg = ipcbeg
      ieend = ipcend
      jebeg = jdcbeg
      jeend = jpcbeg - 1
      kebeg = kdcbeg
      keend = kpcbeg - 1
    CASE (10)
      iebeg = ipcbeg
      ieend = ipcend
      jebeg = jdcbeg
      jeend = jpcbeg - 1
      kebeg = kpcend + 1
      keend = kdcend
    CASE (11)
      iebeg = ipcbeg
      ieend = ipcend
      jebeg = jpcend + 1
      jeend = jdcend
      kebeg = kpcend + 1
      keend = kdcend
    CASE (12)
      iebeg = ipcbeg
      ieend = ipcend
      jebeg = jpcend + 1
      jeend = jdcend
      kebeg = kdcbeg
      keend = kpcbeg - 1
  END SELECT

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetEdgeCellsIndices

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetCornerCellsIndices( region,iLev,icorner, &
                                       icbeg,icend,jcbeg,jcend,kcbeg,kcend )

  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, icorner, icbeg, icend, jcbeg, jcend, kcbeg, kcend
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetCornerCellsIndices',&
  'RFLO_GetDimensions.F90' )

  IF (icorner<0 .OR. icorner>8) &
    CALL ErrorStop( region%global,ERR_VOLUME_CORNERS,&
    __LINE__ )

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )

  SELECT CASE (icorner)
    CASE (1)
      icbeg = idcbeg
      icend = ipcbeg - 1
      jcbeg = jdcbeg
      jcend = jpcbeg - 1
      kcbeg = kdcbeg
      kcend = kpcbeg - 1
    CASE (2)
      icbeg = idcbeg
      icend = ipcbeg - 1
      jcbeg = jdcbeg
      jcend = jpcbeg - 1
      kcbeg = kpcend + 1
      kcend = kdcend
    CASE (3)
      icbeg = idcbeg
      icend = ipcbeg - 1
      jcbeg = jpcend + 1
      jcend = jdcend
      kcbeg = kpcend + 1
      kcend = kdcend
    CASE (4)
      icbeg = idcbeg
      icend = ipcbeg - 1
      jcbeg = jpcend + 1
      jcend = jdcend
      kcbeg = kdcbeg
      kcend = kpcbeg - 1
    CASE (5)
      icbeg = ipcend + 1
      icend = idcend
      jcbeg = jdcbeg
      jcend = jpcbeg - 1
      kcbeg = kdcbeg
      kcend = kpcbeg - 1
    CASE (6)
      icbeg = ipcend + 1
      icend = idcend
      jcbeg = jdcbeg
      jcend = jpcbeg - 1
      kcbeg = kpcend + 1
      kcend = kdcend
    CASE (7)
      icbeg = ipcend + 1
      icend = idcend
      jcbeg = jpcend + 1
      jcend = jdcend
      kcbeg = kpcend + 1
      kcend = kdcend
    CASE (8)
      icbeg = ipcend + 1
      icend = idcend
      jcbeg = jpcend + 1
      jcend = jdcend
      kcbeg = kdcbeg
      kcend = kpcbeg - 1
  END SELECT

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetCornerCellsIndices

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetDimensions.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.11  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/02/03 19:20:46  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.7  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.4  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.3  2002/01/02 16:20:18  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.2  2001/12/19 23:09:20  jblazek
! Added routines to read grid and solution.
!
! Revision 1.1  2001/12/11 21:59:28  jblazek
! memory allocation added.
!
!******************************************************************************














