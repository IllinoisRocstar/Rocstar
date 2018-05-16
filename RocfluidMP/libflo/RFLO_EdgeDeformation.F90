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
! Purpose: calculate node displacements on those edges whose end points have
!          moved, but the associated boundaries were not updated yet (finest
!          grid only).
!
! Description: points along an edge are shifted using 1-D linear transfinite
!              interpolation (TFI).
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOld     = grid from previous time step.
!
! Output: edgeMoved = flag if discretization at an edge was changed
!         dNode     = updated deformations at edges.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************
!
! $Id: RFLO_EdgeDeformation.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_EdgeDeformation( region,boundMoved,edgeMoved, &
                                 arcLen12,arcLen34,arcLen56,xyzOld,dNode )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint1d
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iEdge, ind

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, l1c, l2c
  INTEGER :: indBeg, indEnd, ijkN, ijkN1, ijkNBeg, ijkNEnd, iNOff, ijNOff
  INTEGER :: switch(12,9)

  REAL(RFREAL) :: arcLen, ds, s, dN(3), dNBeg(3), dNEnd(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_EdgeDeformation',&
  'RFLO_EdgeDeformation.F90' )

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set edge switch -------------------------------------------------------------
! switch(:,1) = begins at boundary
! switch(:,2) = ends on boundary
! switch(:,3) = right boundary
! switch(:,4) = left boundary
! switch(:,5) = direction (from-to boundary)
! switch(:,6) = start index
! switch(:,7) = end index
! switch(:,8) = constant index in 1st direction
! switch(:,9) = constant index in 2nd direction

  switch( 1,:) = (/5, 6, 1, 3, 56, kpnbeg, kpnend, ipnbeg, jpnbeg/)
  switch( 2,:) = (/3, 4, 1, 6, 34, jpnbeg, jpnend, kpnend, ipnbeg/)
  switch( 3,:) = (/5, 6, 1, 4, 56, kpnbeg, kpnend, ipnbeg, jpnend/)
  switch( 4,:) = (/3, 4, 1, 5, 34, jpnbeg, jpnend, kpnbeg, ipnbeg/)
  switch( 5,:) = (/5, 6, 2, 3, 56, kpnbeg, kpnend, ipnend, jpnbeg/)
  switch( 6,:) = (/3, 4, 2, 6, 34, jpnbeg, jpnend, kpnend, ipnend/)
  switch( 7,:) = (/5, 6, 2, 4, 56, kpnbeg, kpnend, ipnend, jpnend/)
  switch( 8,:) = (/3, 4, 2, 5, 34, jpnbeg, jpnend, kpnbeg, ipnend/)
  switch( 9,:) = (/1, 2, 3, 5, 12, ipnbeg, ipnend, jpnbeg, kpnbeg/)
  switch(10,:) = (/1, 2, 3, 6, 12, ipnbeg, ipnend, jpnbeg, kpnend/)
  switch(11,:) = (/1, 2, 4, 5, 12, ipnbeg, ipnend, jpnend, kpnbeg/)
  switch(12,:) = (/1, 2, 4, 6, 12, ipnbeg, ipnend, jpnend, kpnend/)

! edge movement flag ----------------------------------------------------------

  edgeMoved(:) = .false.

  IF (boundMoved(1)) THEN
    edgeMoved( 1) = .true.; edgeMoved( 2) = .true.
    edgeMoved( 3) = .true.; edgeMoved( 4) = .true.
  ENDIF
  IF (boundMoved(2)) THEN
    edgeMoved( 5) = .true.; edgeMoved( 6) = .true.
    edgeMoved( 7) = .true.; edgeMoved( 8) = .true.
  ENDIF
  IF (boundMoved(3)) THEN
    edgeMoved( 1) = .true.; edgeMoved( 5) = .true.
    edgeMoved( 9) = .true.; edgeMoved(10) = .true.
  ENDIF
  IF (boundMoved(4)) THEN
    edgeMoved( 3) = .true.; edgeMoved( 7) = .true.
    edgeMoved(11) = .true.; edgeMoved(12) = .true.
  ENDIF
  IF (boundMoved(5)) THEN
    edgeMoved( 4) = .true.; edgeMoved( 8) = .true.
    edgeMoved( 9) = .true.; edgeMoved(11) = .true.
  ENDIF
  IF (boundMoved(6)) THEN
    edgeMoved( 2) = .true.; edgeMoved( 6) = .true.
    edgeMoved(10) = .true.; edgeMoved(12) = .true.
  ENDIF

! loop over all 12 edges ------------------------------------------------------

  DO iEdge=1,12
    IF ((boundMoved(switch(iEdge,1)) .OR. boundMoved(switch(iEdge,2))) .AND. &
        ((.NOT.boundMoved(switch(iEdge,3))) .OR. &
         (.NOT.boundMoved(switch(iEdge,4)))) .AND. &
         (.NOT.edgeMoved(iEdge))) THEN

      edgeMoved(iEdge) = .true.

      ds     = 0._RFREAL
      indBeg = switch(iEdge,6)
      indEnd = switch(iEdge,7)
      l1c    = switch(iEdge,8)
      l2c    = switch(iEdge,9)
      DO ind=indBeg+1,indEnd-1
        IF (switch(iEdge,5) == 12) THEN
          ijkN     = IndIJK(ind   ,l1c,l2c,iNOff,ijNOff)
          ijkN1    = IndIJK(ind-1 ,l1c,l2c,iNOff,ijNOff)
          ijkNBeg  = IndIJK(indBeg,l1c,l2c,iNOff,ijNOff)
          ijkNEnd  = IndIJK(indEnd,l1c,l2c,iNOff,ijNOff)
          arcLen   = arcLen12(l1c,l2c)
          dNBeg(:) = dNode(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd)
        ELSE IF (switch(iEdge,5) == 34) THEN
          ijkN     = IndIJK(l2c,ind   ,l1c,iNOff,ijNOff)
          ijkN1    = IndIJK(l2c,ind-1 ,l1c,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l2c,indBeg,l1c,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l2c,indEnd,l1c,iNOff,ijNOff)
          arcLen   = arcLen34(l1c,l2c)
          dNBeg(:) = dNode(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd)
        ELSE IF (switch(iEdge,5) == 56) THEN
          ijkN     = IndIJK(l1c,l2c,ind   ,iNOff,ijNOff)
          ijkN1    = IndIJK(l1c,l2c,ind-1 ,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l1c,l2c,indBeg,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l1c,l2c,indEnd,iNOff,ijNOff)
          arcLen   = arcLen56(l1c,l2c)
          dNBeg(:) = dNode(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd)
        ENDIF
        ds = ds + SQRT((xyzOld(XCOORD,ijkN)-xyzOld(XCOORD,ijkN1))**2 + &
                       (xyzOld(YCOORD,ijkN)-xyzOld(YCOORD,ijkN1))**2 + &
                       (xyzOld(ZCOORD,ijkN)-xyzOld(ZCOORD,ijkN1))**2)
        s  = ds/arcLen

        CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
        dNode(:,ijkN) = dN(:)
      ENDDO   ! i

    ENDIF     ! boundMoved
  ENDDO       ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_EdgeDeformation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_EdgeDeformation.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/05/27 01:53:41  wasistho
! added rflo_gridremesh
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
! Revision 1.3  2003/03/14 22:05:10  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
!******************************************************************************







