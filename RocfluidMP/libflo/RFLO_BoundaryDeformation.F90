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
! Purpose: calculate node displacements on those boundaries whose edges
!          have moved but which were not marked as moving (finest grid only).
!
! Description: none.
!
! Input: region     = grid dimensions
!        boundMoved = flag for boundaries of a region which have moved
!        edgeMoved  = flag for edges whose nodes have moved
!        arcLen12   = arclength between i=const. boundaries for each j, k
!        arcLen34   = arclength between j=const. boundaries for each k, i
!        arcLen56   = arclength between k=const. boundaries for each i, j
!        xyzOld     = grid from previous time step.
!
! Output: dNode = updated deformations at boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************
!
! $Id: RFLO_BoundaryDeformation.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_BoundaryDeformation( region,boundMoved,edgeMoved, &
                                     arcLen12,arcLen34,arcLen56,  &
                                     xyzOld,dNode )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes, &
                            RFLO_Tfint2d
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
  INTEGER :: iBound, l1, l2

! ... local variables
  INTEGER :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: l1b, l1e, l2b, l2e, lc, ijkN, ijkE(4), ijkEm(4), iNOff, ijNOff
  INTEGER :: switch(6,9)

  LOGICAL :: sum12

  REAL(RFREAL) :: arcLen(4), ds(4), s(4)
  REAL(RFREAL) :: corner(3,8), e1(3), e2(3), e3(3), e4(3), &
                  p1(3), p2(3), p3(3), p4(3), dN(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BoundaryDeformation',&
       'RFLO_BoundaryDeformation.F90' )

! get dimensions --------------------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! set boundary switch ---------------------------------------------------------
! switch(:,1-4) = numbers of the 4 edges of a boundary
! switch(:,5-6) = first/last index in l1-direction
! switch(:,7-8) = first/last index in l2-direction
! switch(:,  9) = constant index

  switch(1,:) = (/ 1,  2,  3,  4, jpnbeg, jpnend, kpnbeg, kpnend, ipnbeg/)
  switch(2,:) = (/ 5,  6,  7,  8, jpnbeg, jpnend, kpnbeg, kpnend, ipnend/)
  switch(3,:) = (/ 1,  5,  9, 10, kpnbeg, kpnend, ipnbeg, ipnend, jpnbeg/)
  switch(4,:) = (/ 3,  7, 11, 12, kpnbeg, kpnend, ipnbeg, ipnend, jpnend/)
  switch(5,:) = (/ 4,  8,  9, 11, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg/)
  switch(6,:) = (/ 2,  6, 10, 12, ipnbeg, ipnend, jpnbeg, jpnend, kpnend/)

! store displacements at corners ----------------------------------------------

  corner(:,1) = dNode(:,IndIJK(ipnbeg,jpnbeg,kpnbeg,iNOff,ijNOff))
  corner(:,2) = dNode(:,IndIJK(ipnbeg,jpnbeg,kpnend,iNOff,ijNOff))
  corner(:,3) = dNode(:,IndIJK(ipnbeg,jpnend,kpnend,iNOff,ijNOff))
  corner(:,4) = dNode(:,IndIJK(ipnbeg,jpnend,kpnbeg,iNOff,ijNOff))
  corner(:,5) = dNode(:,IndIJK(ipnend,jpnbeg,kpnbeg,iNOff,ijNOff))
  corner(:,6) = dNode(:,IndIJK(ipnend,jpnbeg,kpnend,iNOff,ijNOff))
  corner(:,7) = dNode(:,IndIJK(ipnend,jpnend,kpnend,iNOff,ijNOff))
  corner(:,8) = dNode(:,IndIJK(ipnend,jpnend,kpnbeg,iNOff,ijNOff))

! move nodes on boundaries with active edges ----------------------------------

  DO iBound=1,6
    IF ((.NOT.boundMoved(iBound)) .AND. &
        (edgeMoved(switch(iBound,1)) .OR. edgeMoved(switch(iBound,2)) .OR. &
         edgeMoved(switch(iBound,3)) .OR. edgeMoved(switch(iBound,4)))) THEN

      l1b = switch(iBound,5)
      l1e = switch(iBound,6)
      l2b = switch(iBound,7)
      l2e = switch(iBound,8)
      lc  = switch(iBound,9)

      IF (iBound == 1) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,4)
        p3(:) = corner(:,3)
        p4(:) = corner(:,2)
      ELSE IF (iBound == 2) THEN
        p1(:) = corner(:,5)
        p2(:) = corner(:,8)
        p3(:) = corner(:,7)
        p4(:) = corner(:,6)
      ELSE IF (iBound == 3) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,2)
        p3(:) = corner(:,6)
        p4(:) = corner(:,5)
      ELSE IF (iBound == 4) THEN
        p1(:) = corner(:,4)
        p2(:) = corner(:,3)
        p3(:) = corner(:,7)
        p4(:) = corner(:,8)
      ELSE IF (iBound == 5) THEN
        p1(:) = corner(:,1)
        p2(:) = corner(:,5)
        p3(:) = corner(:,8)
        p4(:) = corner(:,4)
      ELSE IF (iBound == 6) THEN
        p1(:) = corner(:,2)
        p2(:) = corner(:,6)
        p3(:) = corner(:,7)
        p4(:) = corner(:,3)
      ENDIF

      ds(1:2) = 0._RFREAL
      DO l2=l2b+1,l2e-1

        sum12   = .true.
        ds(3:4) = 0._RFREAL
        DO l1=l1b+1,l1e-1
          IF (iBound==1 .OR. iBound==2) THEN
            ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(lc,jpnbeg,l2    ,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(lc,jpnbeg,l2-1  ,iNOff,ijNOff)
            ijkE(2)   = IndIJK(lc,jpnend,l2    ,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(lc,jpnend,l2-1  ,iNOff,ijNOff)
            ijkE(3)   = IndIJK(lc,l1    ,kpnbeg,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(lc,l1-1  ,kpnbeg,iNOff,ijNOff)
            ijkE(4)   = IndIJK(lc,l1    ,kpnend,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(lc,l1-1  ,kpnend,iNOff,ijNOff)
            arcLen(1) = arcLen56(lc,jpnbeg)
            arcLen(2) = arcLen56(lc,jpnend)
            arcLen(3) = arcLen34(kpnbeg,lc)
            arcLen(4) = arcLen34(kpnend,lc)
          ELSE IF (iBound==3 .OR. iBound==4) THEN
            ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
            ijkE(1)   = IndIJK(l2    ,lc,kpnbeg,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(l2-1  ,lc,kpnbeg,iNOff,ijNOff)
            ijkE(2)   = IndIJK(l2    ,lc,kpnend,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(l2-1  ,lc,kpnend,iNOff,ijNOff)
            ijkE(3)   = IndIJK(ipnbeg,lc,l1    ,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(ipnbeg,lc,l1-1  ,iNOff,ijNOff)
            ijkE(4)   = IndIJK(ipnend,lc,l1    ,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(ipnend,lc,l1-1  ,iNOff,ijNOff)
            arcLen(1) = arclen12(lc,kpnbeg)
            arcLen(2) = arclen12(lc,kpnend)
            arcLen(3) = arcLen56(ipnbeg,lc)
            arcLen(4) = arcLen56(ipnend,lc)
          ELSE IF (iBound==5 .OR. iBound==6) THEN
            ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
            ijkE(1)   = IndIJK(ipnbeg,l2    ,lc,iNOff,ijNOff)
            ijkEm(1)  = IndIJK(ipnbeg,l2-1  ,lc,iNOff,ijNOff)
            ijkE(2)   = IndIJK(ipnend,l2    ,lc,iNOff,ijNOff)
            ijkEm(2)  = IndIJK(ipnend,l2-1  ,lc,iNOff,ijNOff)
            ijkE(3)   = IndIJK(l1    ,jpnbeg,lc,iNOff,ijNOff)
            ijkEm(3)  = IndIJK(l1-1  ,jpnbeg,lc,iNOff,ijNOff)
            ijkE(4)   = IndIJK(l1    ,jpnend,lc,iNOff,ijNOff)
            ijkEm(4)  = IndIJK(l1-1  ,jpnend,lc,iNOff,ijNOff)
            arcLen(1) = arcLen34(lc,ipnbeg)
            arcLen(2) = arcLen34(lc,ipnend)
            arcLen(3) = arcLen12(jpnbeg,lc)
            arcLen(4) = arcLen12(jpnend,lc)
          ENDIF
          IF (sum12) THEN
            ds(1) = ds(1) + &
                    SQRT((xyzOld(XCOORD,ijkE(1))-xyzOld(XCOORD,ijkEm(1)))**2 + &
                         (xyzOld(YCOORD,ijkE(1))-xyzOld(YCOORD,ijkEm(1)))**2 + &
                         (xyzOld(ZCOORD,ijkE(1))-xyzOld(ZCOORD,ijkEm(1)))**2)
            ds(2) = ds(2) + &
                    SQRT((xyzOld(XCOORD,ijkE(2))-xyzOld(XCOORD,ijkEm(2)))**2 + &
                         (xyzOld(YCOORD,ijkE(2))-xyzOld(YCOORD,ijkEm(2)))**2 + &
                         (xyzOld(ZCOORD,ijkE(2))-xyzOld(ZCOORD,ijkEm(2)))**2)
            sum12 = .false.
          ENDIF
          ds(3) = ds(3) + &
                  SQRT((xyzOld(XCOORD,ijkE(3))-xyzOld(XCOORD,ijkEm(3)))**2 + &
                       (xyzOld(YCOORD,ijkE(3))-xyzOld(YCOORD,ijkEm(3)))**2 + &
                       (xyzOld(ZCOORD,ijkE(3))-xyzOld(ZCOORD,ijkEm(3)))**2)
          ds(4) = ds(4) + &
                  SQRT((xyzOld(XCOORD,ijkE(4))-xyzOld(XCOORD,ijkEm(4)))**2 + &
                       (xyzOld(YCOORD,ijkE(4))-xyzOld(YCOORD,ijkEm(4)))**2 + &
                       (xyzOld(ZCOORD,ijkE(4))-xyzOld(ZCOORD,ijkEm(4)))**2)
          s(:)  = ds(:)/arcLen(:)
          e1(:) = dNode(:,ijkE(1))
          e2(:) = dNode(:,ijkE(2))
          e3(:) = dNode(:,ijkE(3))
          e4(:) = dNode(:,ijkE(4))
          CALL RFLO_Tfint2d( s(1),s(2),s(3),s(4),e1,e2,e3,e4,p1,p2,p3,p4,dN )
          dNode(:,ijkN) = dN(:)
        ENDDO  ! l1
      ENDDO    ! l2

    ENDIF      ! .NOT.boundMoved & edgeMoved
  ENDDO        ! iBound

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BoundaryDeformation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_BoundaryDeformation.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:15  wasistho
! lower to upper case
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.4  2003/05/06 20:05:38  jblazek
! Corrected bug in grid motion (corner "averaging").
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







