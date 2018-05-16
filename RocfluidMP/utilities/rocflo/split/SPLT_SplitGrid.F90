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
! Purpose: split single grid into multiple regions.
!
! Description: none.
!
! Input: splitDirection = direction (i,j,k) in which the grid is to be splitted
!        regionsOld     = single grid and topology.
!
! Output: regionsNew = splitted grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SPLT_SplitGrid.F90,v 1.4 2008/12/06 08:44:51 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SplitGrid( splitDirection,regionsOld,regionsNew )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE SPLT_ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: splitDirection

  TYPE(t_region), POINTER :: regionsOld(:), regionsNew(:)

! ... loop variables
  INTEGER :: iReg, i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iNOff, ijNOff
  INTEGER :: idnbegOld, idnendOld, jdnbegOld, jdnendOld, kdnbegOld, kdnendOld, &
             iNOffOld, ijNOffOld, errorFlag
  INTEGER :: ibn, ien, ijkN, ijkNOld, ii, jj, kk, nCells, nPatches

  TYPE(t_global), POINTER :: globalNew
  TYPE(t_grid), POINTER   :: gridOld, gridNew

!******************************************************************************

  globalNew => regionsNew(1)%global

  CALL RegisterFunction( globalNew,'SplitGrid',&
  'SPLT_SplitGrid.F90' )

! allocate memory for new grid; set dimensions --------------------------------

  gridOld  => regionsOld(1)%levels(1)%grid
  nPatches =  regionsOld(1)%nPatches + 2     ! max. number of patches

  DO iReg=1,globalNew%nRegions

    ALLOCATE( regionsNew(iReg)%levels(1),stat=errorFlag )
    ALLOCATE( regionsNew(iReg)%levels(1)%patches(nPatches),stat=errorFlag )
    globalNew%error = errorFlag
    IF (globalNew%error /= 0) CALL ErrorStop( globalNew,ERR_ALLOCATE,__LINE__ )

    gridNew => regionsNew(iReg)%levels(1)%grid

    IF (splitDirection == 1) THEN
      gridNew%ipc = gridOld%ipc/globalNew%nRegions
      gridNew%jpc = gridOld%jpc
      gridNew%kpc = gridOld%kpc
      IF (iReg == globalNew%nRegions) THEN
        gridNew%ipc = gridOld%ipc - gridNew%ipc*(globalNew%nRegions-1)
      ENDIF
    ELSE IF (splitDirection == 2) THEN
      gridNew%ipc = gridOld%ipc
      gridNew%jpc = gridOld%jpc/globalNew%nRegions
      gridNew%kpc = gridOld%kpc
      IF (iReg == globalNew%nRegions) THEN
        gridNew%jpc = gridOld%jpc - gridNew%jpc*(globalNew%nRegions-1)
      ENDIF
    ELSE
      gridNew%ipc = gridOld%ipc
      gridNew%jpc = gridOld%jpc
      gridNew%kpc = gridOld%kpc/globalNew%nRegions
      IF (iReg == globalNew%nRegions) THEN
        gridNew%kpc = gridOld%kpc - gridNew%kpc*(globalNew%nRegions-1)
      ENDIF
    ENDIF
    regionsNew(iReg)%nGridLevels = regionsOld(1)%nGridLevels
    regionsNew(iReg)%nDumCells   = 0
    regionsNew(iReg)%nPatches    = 0

    CALL RFLO_GetDimensDummyNodes( regionsNew(iReg),1,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regionsNew(iReg),1,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
    ALLOCATE( gridNew%xyz(3,ibn:ien),stat=errorFlag )
    globalNew%error = errorFlag
    IF (globalNew%error /= 0) CALL ErrorStop( globalNew,ERR_ALLOCATE,__LINE__ )

  ENDDO   ! iReg

! copy grid to new regions ----------------------------------------------------

  gridOld => regionsOld(1)%levels(1)%grid
  CALL RFLO_GetDimensDummyNodes( regionsOld(1),1,idnbegOld,idnendOld, &
                                 jdnbegOld,jdnendOld,kdnbegOld,kdnendOld )
  CALL RFLO_GetNodeOffset( regionsOld(1),1,iNOffOld,ijNOffOld )

  DO iReg=1,globalNew%nRegions

    gridNew => regionsNew(iReg)%levels(1)%grid
    CALL RFLO_GetDimensDummyNodes( regionsNew(iReg),1,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( regionsNew(iReg),1,iNOff,ijNOff )

    DO k=kdnbeg,kdnend
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
          IF (splitDirection == 1) THEN
            nCells  = gridOld%ipc/globalNew%nRegions
            ii      = i + (iReg-1)*nCells
            ijkNOld = IndIJK(ii,j,k,iNOffOld,ijNOffOld)
          ELSE IF (splitDirection == 2) THEN
            nCells  = gridOld%jpc/globalNew%nRegions
            jj      = j + (iReg-1)*nCells
            ijkNOld = IndIJK(i,jj,k,iNOffOld,ijNOffOld)
          ELSE
            nCells  = gridOld%kpc/globalNew%nRegions
            kk      = k + (iReg-1)*nCells
            ijkNOld = IndIJK(i,j,kk,iNOffOld,ijNOffOld)
          ENDIF
          gridNew%xyz(XCOORD,ijkN) = gridOld%xyz(XCOORD,ijkNOld)
          gridNew%xyz(YCOORD,ijkN) = gridOld%xyz(YCOORD,ijkNOld)
          gridNew%xyz(ZCOORD,ijkN) = gridOld%xyz(ZCOORD,ijkNOld)
        ENDDO
      ENDDO
    ENDDO

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( globalNew )

END SUBROUTINE SplitGrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPLT_SplitGrid.F90,v $
! Revision 1.4  2008/12/06 08:44:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:33:31  wasistho
! rflo_modinterfacessplit to splt_modinterfaces
!
! Revision 1.1  2004/12/03 02:41:15  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:45:30  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.5  2003/03/20 22:32:37  haselbac
! Renamed ModInterfaces
!
! Revision 1.4  2003/03/20 19:45:54  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.3  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/12 21:50:08  jblazek
! Added tool to split single grid into multiple regions.
!
!******************************************************************************







