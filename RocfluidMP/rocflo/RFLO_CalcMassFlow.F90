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
! Purpose: calculate amount of inflow and outflow for a region.
!
! Description: none.
!
! Input: region%levels%mixt        = flow variables
!        region%levels%grid%si/j/k = face vectors (at boundaries)
!
! Output: global%massIn  = inflow mass
!         global%massOut = outflow mass
!
! Notes: mass flow is computed as a sum over all inflow/outflow/injection
!        boundaries of the region.
!
!******************************************************************************
!
! $Id: RFLO_CalcMassFlow.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcMassFlow( region )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            RFLO_GetPatchIndices, RFLO_GetPatchDirection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: bcType, iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkD, ijkN
  INTEGER :: inode, jnode, knode

  LOGICAL :: inflow, outflow

  REAL(RFREAL)          :: sgn, massIn, massOut, rhoua, rhova, rhowa, mass
  REAL(RFREAL), POINTER :: cv(:,:), sFace(:,:)

  TYPE(t_patch), POINTER :: patch

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcMassFlow',&
  'RFLO_CalcMassFlow.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev =  region%currLevel

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv => region%levels(iLev)%mixt%cv

  massIn  = region%global%massIn
  massOut = region%global%massOut

! loop over patches -----------------------------------------------------------

  DO iPatch=1,region%nPatches

    patch  => region%levels(iLev)%patches(iPatch)
    bcType =  patch%bcType
    lbound =  patch%lbound

    IF ((bcType>=BC_INFLOW    .AND. bcType<=BC_INFLOW   +BC_RANGE) .OR. &
        (bcType>=BC_OUTFLOW   .AND. bcType<=BC_OUTFLOW  +BC_RANGE) .OR. &
        (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) .OR. &
        (bcType>=BC_FARFIELD  .AND. bcType<=BC_FARFIELD +BC_RANGE)) THEN

! --- get dimensions of the patch

      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

! --- inflow or outflow?

      inflow  = .false.
      outflow = .false.
      IF (bcType>=BC_INFLOW    .AND. bcType<=BC_INFLOW   +BC_RANGE) inflow  = .true.
      IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) inflow  = .true.
      IF (bcType>=BC_OUTFLOW   .AND. bcType<=BC_OUTFLOW  +BC_RANGE) outflow = .true.
      IF (bcType>=BC_FARFIELD  .AND. bcType<=BC_FARFIELD +BC_RANGE) outflow = .true.

! --- to take the right face vector and make it point outwards

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

! --- get the appropriate face vector

      IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
      IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
      IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! --- loop over all cells of the patch

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC  = IndIJK(i,j,k,iCOff,ijCOff)
            ijkD  = IndIJK(i-idir,j-jdir,k-kdir,iCOff,ijCOff)
            ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            rhoua = 0.5_RFREAL*(cv(CV_MIXT_XMOM,ijkC)+cv(CV_MIXT_XMOM,ijkD))
            rhova = 0.5_RFREAL*(cv(CV_MIXT_YMOM,ijkC)+cv(CV_MIXT_YMOM,ijkD))
            rhowa = 0.5_RFREAL*(cv(CV_MIXT_ZMOM,ijkC)+cv(CV_MIXT_ZMOM,ijkD))
            mass  = sgn*sFace(XCOORD,ijkN)*rhoua + &
                    sgn*sFace(YCOORD,ijkN)*rhova + &
                    sgn*sFace(ZCOORD,ijkN)*rhowa
            IF (inflow ) massIn  = massIn  - mass
            IF (outflow) massOut = massOut + mass
          ENDDO
        ENDDO
      ENDDO

    ENDIF   ! inflow/outflow

  ENDDO     ! iPatch

! finalize --------------------------------------------------------------------

  region%global%massIn  = massIn
  region%global%massOut = massOut

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcMassFlow

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcMassFlow.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/20 20:06:41  wasistho
! enable mass flow computation at farfield boundary
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.16  2003/11/20 16:40:37  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.12  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.11  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.10  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.9  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.8  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
! Revision 1.7  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.5  2002/02/01 22:17:38  jblazek
! Change addressing of face vectors at block boundaries.
!
! Revision 1.4  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.3  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.2  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







