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
! Purpose: extrapolate gradients from the interior domain or at the patch 
!          to the dummy points depending on the boundary condition.
!
! Description: dummy gradients of physical boundary copied from interior
!              dummy gradients of symmetry boundary mirrored from interior
!              outerst dummy gradients of connecting bnd copied from interior
!
! Input: region  = data of current region.
!        patch   = current patch.
!        iBegV, iEndV = begin and end var index
!        iBegG, iEndG = begin and end gradient index
!  
! Output: gradi, gradj, gradk = gradients at dummy points of the current patch
!                               plane 
!
! Notes: many parameters are passed as subroutine arguments to minimize
!        repetition of work
!
!******************************************************************************
!
! $Id: RFLO_CalcGradDummy.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradDummy( region,patch,iBegV,iEndV,iBegG,iEndG, &
                               gradi,gradj,gradk )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices,RFLO_GetPatchDirection, &
                            RFLO_CalcGradDummyPhys,RFLO_CalcGradDummyConn, &
                            RFLO_CalcGradDummySymm
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

  INTEGER :: iBegV, iEndV, iBegG, iEndG

  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)

! ... local variables
  INTEGER :: iLev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  INTEGER :: inode, jnode, knode, idir, jdir, kdir

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradDummy',&
       'RFLO_CalcGradDummy.F90' )

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
      (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE)) THEN

! - complete dummy layers by zeroth extrapolation 

    CALL RFLO_CalcGradDummyPhys( region,lbound, &
                                 idir  ,jdir  ,kdir  , &
                                 indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                 iBegV ,iEndV ,iBegG ,iEndG , &
                                 gradi ,gradj ,gradk )

  ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN

! - complete dummy layers by linear extrapolation

    CALL RFLO_CalcGradDummySymm( region,lbound, &
                                  idir  ,jdir  ,kdir  , &
                                  inode ,jnode ,knode , &
                                  ibeg  ,iend  ,jbeg  ,jend  ,kbeg  ,kend  , &
                                  indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                  iBegV ,iEndV ,iBegG ,iEndG , &
                                  gradi ,gradj ,gradk )

  ELSE IF ((bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
           (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) .OR. &
           (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE)) THEN

! - complete the outerst dummy layers by zeroth extrapolation

    CALL RFLO_CalcGradDummyConn( region,lbound, &
                                 idir  ,jdir  ,kdir  , &
                                 indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                 iBegV ,iEndV ,iBegG ,iEndG , &
                                 gradi ,gradj ,gradk )

  ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! - everything else

  ELSE

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcGradDummy

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradDummy.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/09/02 22:58:54  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:43:11  wasistho
! made compliant with CODING RULE
!
!******************************************************************************







