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
! Purpose: Extrapolate face-widths from the interior domain or at the patch 
!          to the dummy faces depending on the boundary condition.
!
! Description: Dummy values of physical and symmetry boundary are mirrored
!              Dummy values of connecting bnd has been computed in prev. step 
!
! Input: region  = data of current region.
!        patch   = current patch.
!  
! Output: Face widths at dummy points stored in workI, workJ, workK.
!
! Notes: Mother routine = TURB_FloFaceWidth.
!
!******************************************************************************
!
! $Id: TURB_floFaceWidthDummy.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloFaceWidthDummy( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection
  USE TURB_ModInterfaces, ONLY : TURB_FloFaceWidthDummyPhys, &
                                 TURB_FloFaceWidthDummyConn
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global


  INTEGER :: iLev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  INTEGER :: inode, jnode, knode, idir, jdir, kdir

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floFaceWidthDummy.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloFaceWidthDummy',&
  'TURB_floFaceWidthDummy.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound
  bcType = patch%bcType

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

! take the right face vector and make it point outwards ----------------------

  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

  indBeg = ibeg + inode
  jndBeg = jbeg + jnode
  kndBeg = kbeg + knode
  indEnd = iend + inode + 1 - ABS(idir)
  jndEnd = jend + jnode + 1 - ABS(jdir)
  kndEnd = kend + knode + 1 - ABS(kdir)
 
! call dummy extrapolation routines depending on bcType -----------------------

  IF ((bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR. &
      (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
      (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE) .OR. &
      (bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
      (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR. &
      (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
      (bcType>=BC_SYMMETRY   .AND. bcType<=BC_SYMMETRY  +BC_RANGE)) THEN

! - complete dummy layers by mirroring

    CALL TURB_FloFaceWidthDummyPhys( region,lbound,idir,jdir,kdir, &
                                  indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )

  ELSE IF ((bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
           (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
           (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE)) THEN

! - complete the right outerst dummy layers by zeroth extrapolation

    IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
      CALL TURB_FloFaceWidthDummyConn( region,lbound,idir,jdir,kdir, &
                                  indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )
    ENDIF

  ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! - everything else

  ELSE

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloFaceWidthDummy

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_floFaceWidthDummy.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:48:19  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







