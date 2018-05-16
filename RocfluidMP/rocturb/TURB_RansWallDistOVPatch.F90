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
! Purpose: compute distance of each cell centers to the nearest no-slip wall
!
! Description: in this routine we search for no-slip wall, for each no-slip
!              wall patch, we loop over regions and compute the minimum
!              wall distance. In this way, global wall distance is obtained,
!              assuming open view (OV) from cell centers to the nearest wall.
!
! Input: region = data of current region
!
! Output: region%levels%turb%lens = turbulence length scale.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansWallDistOVPatch.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansWallDistOVPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : CentroidHexa, &
                            RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, l, m

! ... local variables
  INTEGER :: bcType, ibeg, iend, jbeg, jend, kbeg, kend, corner(8)
  REAL(RFREAL) :: csfx, csfy, csfz
  REAL(RFREAL), POINTER :: xyz(:,:)

#ifdef RFLO
  INTEGER :: iLev, lbound, inode, jnode, knode, ia, ja, ka
  INTEGER :: idir, jdir, kdir, iCOff, ijCOff, iNOff, ijNOff
  REAL(RFREAL) :: sgn, xyzHexa(NDIR,8)
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER :: fc(:,:)
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_RansWallDistOVPatch',&
  'TURB_RansWallDistOVPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  bcType =  patch%bcType

#ifdef RFLO
  iLev   = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  lbound =  patch%lbound
  xyz    => region%levels(iLev)%grid%xyz

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector and grid speed

  IF (lbound==1 .OR. lbound==2) THEN
    ia=0
    ja=1
    ka=1
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    ia=1
    ja=0
    ka=1
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    ia=1
    ja=1
    ka=0
  ENDIF
#endif
#ifdef RFLU
  xyz  => region%grid%xyz
  fc   => region%grid%fc
  ibeg =  1
  iend =  patch%nBFaces
#endif

  IF ((bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
      (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE)) THEN

#ifdef RFLO
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          corner(1) = IndIJK(i+inode   ,j+jnode   ,k+knode   ,iNOff,ijNOff)
          corner(2) = IndIJK(i+inode   ,j+jnode   ,k+knode+ka,iNOff,ijNOff)
          corner(3) = IndIJK(i+inode   ,j+jnode+ja,k+knode+ka,iNOff,ijNOff)
          corner(4) = IndIJK(i+inode   ,j+jnode+ja,k+knode   ,iNOff,ijNOff)
          corner(5) = IndIJK(i+inode+ia,j+jnode   ,k+knode   ,iNOff,ijNOff)
          corner(6) = IndIJK(i+inode+ia,j+jnode   ,k+knode+ka,iNOff,ijNOff)
          corner(7) = IndIJK(i+inode+ia,j+jnode+ja,k+knode+ka,iNOff,ijNOff)
          corner(8) = IndIJK(i+inode+ia,j+jnode+ja,k+knode   ,iNOff,ijNOff)

          DO m = 1,8
            DO l = 1,3
              xyzHexa(l,m) = xyz(l,corner(m))
            ENDDO
          ENDDO

          CALL CentroidHexa( xyzHexa,csfx,csfy,csfz )
#endif
#ifdef RFLU
    DO i=ibeg,iend

          csfx = fc(XCOORD,i)
          csfy = fc(YCOORD,i)
          csfz = fc(ZCOORD,i)
#endif
          region%global%turbWorkDim = region%global%turbWorkDim + 1
          region%global%turbWork1D(region%global%turbWorkDim) = csfx
          region%global%turbWorkDim = region%global%turbWorkDim + 1
          region%global%turbWork1D(region%global%turbWorkDim) = csfy
          region%global%turbWorkDim = region%global%turbWorkDim + 1
          region%global%turbWork1D(region%global%turbWorkDim) = csfz
#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k of patch coordinate
#endif
#ifdef RFLU
    ENDDO       ! i
#endif
  ENDIF         ! bcType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_RansWallDistOVPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansWallDistOVPatch.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/23 03:35:00  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.4  2004/01/22 04:02:23  wasistho
! add wall-injection bctype for wall distance
!
! Revision 1.3  2003/10/17 20:21:59  wasistho
! do-loop i.o. dimension range
!
! Revision 1.2  2003/10/07 20:32:24  wasistho
! turbWork2D to turbWork1D
!
!
!******************************************************************************







