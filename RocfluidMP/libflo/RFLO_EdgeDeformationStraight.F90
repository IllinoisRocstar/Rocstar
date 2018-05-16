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
! Purpose: calculate node displacements on straight edges whose points have
!          been moved by RFLO_EdgeDeformation (finest grid only).
!
! Description: points along an edge are shifted using 1-D linear transfinite
!              interpolation (TFI).
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        boundFlat  = flag for boundaries of a region which are flat
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOrig    = reference grid which TFI is based on
!        xyzOld     = grid from previous time step.
!
! Output: edgeMoved = flag if discretization at an edge was changed
!         dNode     = updated deformation at edges.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************
!
! $Id: RFLO_EdgeDeformationStraight.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_EdgeDeformationStraight( region,boundMoved,edgeStraight, &
                                   edgeMoved,arcLen12,arcLen34,arcLen56, &
                                   xyzOrig,xyzOld,dNode )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint1d
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6), edgeStraight(12), edgeMoved(12)

  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:), xyzOrig(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iEdge, ind

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, l1c, l2c
  INTEGER :: indBeg, indEnd, ijkN, ijkN1, ijkNBeg, ijkNEnd, iNOff, ijNOff
  INTEGER :: switch(12,9)

  REAL(RFREAL) :: arcLen, ds, s, dN(3), dNBeg(3), dNEnd(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_EdgeDeformationStraight',&
  'RFLO_EdgeDeformationStraight.F90')

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

! loop over all 12 edges ------------------------------------------------------

  DO iEdge=1,12

    IF ((.NOT.boundMoved(switch(iEdge,3))) .AND. &
        (.NOT.boundMoved(switch(iEdge,4))) .AND. &
        (edgeStraight(iEdge) .EQV. .TRUE.) .AND. &
        (edgeMoved(   iEdge) .EQV. .TRUE.)) THEN

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
          dNBeg(:) = dNode(:,ijkNBeg) + xyzOld(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd) + xyzOld(:,ijkNEnd)
        ELSE IF (switch(iEdge,5) == 34) THEN
          ijkN     = IndIJK(l2c,ind   ,l1c,iNOff,ijNOff)
          ijkN1    = IndIJK(l2c,ind-1 ,l1c,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l2c,indBeg,l1c,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l2c,indEnd,l1c,iNOff,ijNOff)
          arcLen   = arcLen34(l1c,l2c)
          dNBeg(:) = dNode(:,ijkNBeg) + xyzOld(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd) + xyzOld(:,ijkNEnd)
        ELSE IF (switch(iEdge,5) == 56) THEN
          ijkN     = IndIJK(l1c,l2c,ind   ,iNOff,ijNOff)
          ijkN1    = IndIJK(l1c,l2c,ind-1 ,iNOff,ijNOff)
          ijkNBeg  = IndIJK(l1c,l2c,indBeg,iNOff,ijNOff)
          ijkNEnd  = IndIJK(l1c,l2c,indEnd,iNOff,ijNOff)
          arcLen   = arcLen56(l1c,l2c)
          dNBeg(:) = dNode(:,ijkNBeg) + xyzOld(:,ijkNBeg)
          dNEnd(:) = dNode(:,ijkNEnd) + xyzOld(:,ijkNEnd)
        ENDIF
        ds = ds + SQRT((xyzOrig(XCOORD,ijkN)-xyzOrig(XCOORD,ijkN1))**2 + &
                       (xyzOrig(YCOORD,ijkN)-xyzOrig(YCOORD,ijkN1))**2 + &
                       (xyzOrig(ZCOORD,ijkN)-xyzOrig(ZCOORD,ijkN1))**2)
        s  = ds/arcLen

        CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
        dNode(:,ijkN) = dN(:) - xyzOld(:,ijkN)
      ENDDO   ! i

    ENDIF     ! boundMoved
  ENDDO       ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_EdgeDeformationStraight

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_EdgeDeformationStraight.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/14 04:35:50  wasistho
! improved computation of straight edges
!
! Revision 1.1  2006/03/12 20:42:53  wasistho
! added RFLO_EdgeDeformationStraight
!
!
!******************************************************************************







