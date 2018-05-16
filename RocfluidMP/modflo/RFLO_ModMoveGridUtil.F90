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
! ******************************************************************************
!
! Purpose: Collection of grid-motion utility routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModMoveGridUtil.F90,v 1.5 2008/12/06 08:44:17 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModMoveGridUtil

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_MoveGridQFlatPatch, &
            RFLO_MoveGridCurvedPatch

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModMoveGridUtil.F90,v $ $Revision: 1.5 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


!******************************************************************************
!
! Purpose: move grid on non-external half-flat patch (curved in one direction) 
!          (finest grid only).
!
! Description: none.
!
! Input: region = current region
!        patch  = current patch
!        iPatch = patch number
!
! Output: xyz = updated grid at this patch.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridQFlatPatch( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, &
                            RFLO_Tfint1d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters

  INTEGER :: iPatch
  TYPE(t_region) :: region
  TYPE(t_patch), POINTER :: patch

! ... loop variables
  INTEGER :: l1, l2

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: ijkN, ijkN1, ijkNb, ijkNe, iNOff, ijNOff, errorFlag
  INTEGER :: l1b, l1e, l2b, l2e, lc, dir1, dir2, k1, k2, switch(6,7)

  REAL(RFREAL) :: arcLen, ds1, s
  REAL(RFREAL) :: dNBeg(XCOORD:ZCOORD), dNEnd(XCOORD:ZCOORD), dN(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzRef(:,:)
  REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
  REAL(RFREAL), ALLOCATABLE :: ds2(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_MoveGridQFlatPatch',&
       'RFLO_ModMoveGridUtil.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )

  xyzRef   => region%levels(iLev)%gridOld%xyzOld
  xyz      => region%levels(iLev)%grid%xyz
  arcLen12 => region%levels(iLev)%grid%arcLen12
  arcLen34 => region%levels(iLev)%grid%arcLen34
  arcLen56 => region%levels(iLev)%grid%arcLen56

! set boundary switch ---------------------------------------------------------
! switch(:,1-2) = l1 and l2 patch direction
! switch(:,3-4) = first/last index in l1-direction
! switch(:,5-6) = first/last index in l2-direction
! switch(:,  7) = constant index

  switch(1,:) = (/JCOORD, KCOORD, jbeg, jend, kbeg, kend, ibeg/)
  switch(2,:) = (/JCOORD, KCOORD, jbeg, jend, kbeg, kend, iend/)
  switch(3,:) = (/KCOORD, ICOORD, kbeg, kend, ibeg, iend, jbeg/)
  switch(4,:) = (/KCOORD, ICOORD, kbeg, kend, ibeg, iend, jend/)
  switch(5,:) = (/ICOORD, JCOORD, ibeg, iend, jbeg, jend, kbeg/)
  switch(6,:) = (/ICOORD, JCOORD, ibeg, iend, jbeg, jend, kend/)

  dir1 = switch(lbound,1)
  dir2 = switch(lbound,2)
  l1b  = switch(lbound,3)
  l1e  = switch(lbound,4)
  l2b  = switch(lbound,5)
  l2e  = switch(lbound,6)
  lc   = switch(lbound,7)

  ALLOCATE( ds2(l1b:l1e), stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

! apply 1D TFI along flat direction of Qflat patch ----------------------------

  ds2(:) = 0._RFREAL
  DO l2=l2b+1,l2e-1
    k2 = l2-l2b+1

    ds1 = 0._RFREAL
    DO l1=l1b+1,l1e-1
      k1 = l1-l1b+1

      IF (lbound==1 .OR. lbound==2) THEN
        IF (patch%dirFlat==dir1) THEN
          ijkN  = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
          ijkN1 = IndIJK(lc,l1-1  ,l2    ,iNOff,ijNOff)
          ijkNb = IndIJK(lc,l1b   ,l2    ,iNOff,ijNOff)
          ijkNe = IndIJK(lc,l1e   ,l2    ,iNOff,ijNOff)
          arcLen = patch%arcLen1(k2)
        ELSEIF (patch%dirFlat==dir2) THEN
          ijkN  = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
          ijkN1 = IndIJK(lc,l1    ,l2-1  ,iNOff,ijNOff)
          ijkNb = IndIJK(lc,l1    ,l2b   ,iNOff,ijNOff)
          ijkNe = IndIJK(lc,l1    ,l2e   ,iNOff,ijNOff)
          arcLen = patch%arcLen2(k1)
        ENDIF
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        IF (patch%dirFlat==dir1) THEN
          ijkN  = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
          ijkN1 = IndIJK(l2    ,lc,l1-1  ,iNOff,ijNOff)
          ijkNb = IndIJK(l2    ,lc,l1b   ,iNOff,ijNOff)
          ijkNe = IndIJK(l2    ,lc,l1e   ,iNOff,ijNOff)
          arcLen = patch%arcLen1(k2)
        ELSEIF (patch%dirFlat==dir2) THEN
          ijkN  = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
          ijkN1 = IndIJK(l2-1  ,lc,l1    ,iNOff,ijNOff)
          ijkNb = IndIJK(l2b   ,lc,l1    ,iNOff,ijNOff)
          ijkNe = IndIJK(l2e   ,lc,l1    ,iNOff,ijNOff)
          arcLen = patch%arcLen2(k1)
        ENDIF
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        IF (patch%dirFlat==dir1) THEN
          ijkN  = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
          ijkN1 = IndIJK(l1-1  ,l2    ,lc,iNOff,ijNOff)
          ijkNb = IndIJK(l1b   ,l2    ,lc,iNOff,ijNOff)
          ijkNe = IndIJK(l1e   ,l2    ,lc,iNOff,ijNOff)
          arcLen = patch%arcLen1(k2)
        ELSEIF (patch%dirFlat==dir2) THEN
          ijkN  = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
          ijkN1 = IndIJK(l1    ,l2-1  ,lc,iNOff,ijNOff)
          ijkNb = IndIJK(l1    ,l2b   ,lc,iNOff,ijNOff)
          ijkNe = IndIJK(l1    ,l2e   ,lc,iNOff,ijNOff)
          arcLen = patch%arcLen2(k1)
        ENDIF  ! dirFlat
      ENDIF    ! lbound

      dNBeg(:) = xyz(:,ijkNb)
      dNEnd(:) = xyz(:,ijkNe)

      IF (patch%dirFlat==dir1) THEN
        ds1 = ds1 + &
             SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN1))**2 + &
                  (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN1))**2 + &
                  (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN1))**2)
        s = ds1/arcLen
      ELSEIF (patch%dirFlat==dir2) THEN 
        ds2(l1) = ds2(l1) + &
             SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN1))**2 + &
                  (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN1))**2 + &
                  (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN1))**2)
        s = ds2(l1)/arcLen
      ENDIF

      CALL RFLO_Tfint1d( s,dNBeg,dNEnd,dN )
      xyz(:,ijkN) = dN(:)

    ENDDO ! l1
  ENDDO   ! l2

! deallocate arrays -----------------------------------------------------------

  DEALLOCATE( ds2, stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MoveGridQFlatPatch


!******************************************************************************
!
! Purpose: move grid on non-external arbitrary-curved-patches
!          (finest grid only).
!
! Description: the same method as RFLO_MgFrameBndDeformation is used, but
!              based on total displacement rather than incremental one. 
!
! Input: region = current region
!        patch  = current patch
!        iPatch = patch number
!
! Output: xyz = updated grid at this patch.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_MoveGridCurvedPatch( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes, &
                            RFLO_Tfint2d
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters

  TYPE(t_region)         :: region
  TYPE(t_patch), POINTER :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: l1, l2

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: l1b, l1e, l2b, l2e, lc, ijkN, ijkE(4), ijkEm(4), iNOff, ijNOff
  INTEGER :: h1, h2, switch(6,5)

  LOGICAL :: sum12
  REAL(RFREAL) :: arcLen(4), ds(4), s(4)
  REAL(RFREAL) :: e1(XCOORD:ZCOORD), e2(XCOORD:ZCOORD), &
                  e3(XCOORD:ZCOORD), e4(XCOORD:ZCOORD), &
                  p1(XCOORD:ZCOORD), p2(XCOORD:ZCOORD), &
                  p3(XCOORD:ZCOORD), p4(XCOORD:ZCOORD), dN(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzRef(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_MoveGridCurvedPatch',&
       'RFLO_ModMoveGridUtil.F90' )

! get dimensions and pointers -------------------------------------------------

  lbound = patch%lbound
  iLev   = 1

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )

  xyzRef => region%levels(iLev)%gridOld%xyzOld
  xyz    => region%levels(iLev)%grid%xyz

! set boundary switch ---------------------------------------------------------
! switch(:,1-2) = first/last index in l1-direction
! switch(:,3-4) = first/last index in l2-direction
! switch(:,  5) = constant index

  switch(1,:) = (/jbeg, jend, kbeg, kend, ibeg/)
  switch(2,:) = (/jbeg, jend, kbeg, kend, iend/)
  switch(3,:) = (/kbeg, kend, ibeg, iend, jbeg/)
  switch(4,:) = (/kbeg, kend, ibeg, iend, jend/)
  switch(5,:) = (/ibeg, iend, jbeg, jend, kbeg/)
  switch(6,:) = (/ibeg, iend, jbeg, jend, kend/)

  l1b = switch(lbound,1)
  l1e = switch(lbound,2)
  l2b = switch(lbound,3)
  l2e = switch(lbound,4)
  lc  = switch(lbound,5)
  h1  = l1e-l1b+1
  h2  = l2e-l2b+1

  p1(:) = xyz(:,patch%corns(1)) - xyzRef(:,patch%corns(1)) 
  p2(:) = xyz(:,patch%corns(4)) - xyzRef(:,patch%corns(4))
  p3(:) = xyz(:,patch%corns(3)) - xyzRef(:,patch%corns(3))
  p4(:) = xyz(:,patch%corns(2)) - xyzRef(:,patch%corns(2))

! obtain arclen along patch edges

  arclen(1) = patch%arclen2(1)
  arclen(2) = patch%arclen2(h1)
  arclen(3) = patch%arclen1(1)
  arclen(4) = patch%arclen1(h2)

! conduct TFI on interior patch surface

  ds(1:2) = 0._RFREAL
  DO l2=l2b+1,l2e-1

    sum12   = .true.
    ds(3:4) = 0._RFREAL
    DO l1=l1b+1,l1e-1
      IF (lbound==1 .OR. lbound==2) THEN
        ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
        ijkE(1)   = IndIJK(lc,jbeg  ,l2    ,iNOff,ijNOff)
        ijkEm(1)  = IndIJK(lc,jbeg  ,l2-1  ,iNOff,ijNOff)
        ijkE(2)   = IndIJK(lc,jend  ,l2    ,iNOff,ijNOff)
        ijkEm(2)  = IndIJK(lc,jend  ,l2-1  ,iNOff,ijNOff)
        ijkE(3)   = IndIJK(lc,l1    ,kbeg  ,iNOff,ijNOff)
        ijkEm(3)  = IndIJK(lc,l1-1  ,kbeg  ,iNOff,ijNOff)
        ijkE(4)   = IndIJK(lc,l1    ,kend  ,iNOff,ijNOff)
        ijkEm(4)  = IndIJK(lc,l1-1  ,kend  ,iNOff,ijNOff)
      ELSE IF (lbound==3 .OR. lbound==4) THEN
        ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
        ijkE(1)   = IndIJK(l2    ,lc,kbeg  ,iNOff,ijNOff)
        ijkEm(1)  = IndIJK(l2-1  ,lc,kbeg  ,iNOff,ijNOff)
        ijkE(2)   = IndIJK(l2    ,lc,kend  ,iNOff,ijNOff)
        ijkEm(2)  = IndIJK(l2-1  ,lc,kend  ,iNOff,ijNOff)
        ijkE(3)   = IndIJK(ibeg  ,lc,l1    ,iNOff,ijNOff)
        ijkEm(3)  = IndIJK(ibeg  ,lc,l1-1  ,iNOff,ijNOff)
        ijkE(4)   = IndIJK(iend  ,lc,l1    ,iNOff,ijNOff)
        ijkEm(4)  = IndIJK(iend  ,lc,l1-1  ,iNOff,ijNOff)
      ELSE IF (lbound==5 .OR. lbound==6) THEN
        ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
        ijkE(1)   = IndIJK(ibeg  ,l2    ,lc,iNOff,ijNOff)
        ijkEm(1)  = IndIJK(ibeg  ,l2-1  ,lc,iNOff,ijNOff)
        ijkE(2)   = IndIJK(iend  ,l2    ,lc,iNOff,ijNOff)
        ijkEm(2)  = IndIJK(iend  ,l2-1  ,lc,iNOff,ijNOff)
        ijkE(3)   = IndIJK(l1    ,jbeg  ,lc,iNOff,ijNOff)
        ijkEm(3)  = IndIJK(l1-1  ,jbeg  ,lc,iNOff,ijNOff)
        ijkE(4)   = IndIJK(l1    ,jend  ,lc,iNOff,ijNOff)
        ijkEm(4)  = IndIJK(l1-1  ,jend  ,lc,iNOff,ijNOff)
      ENDIF
      IF (sum12) THEN
        ds(1) = ds(1) + &
                SQRT((xyzRef(XCOORD,ijkE(1))-xyzRef(XCOORD,ijkEm(1)))**2 + &
                     (xyzRef(YCOORD,ijkE(1))-xyzRef(YCOORD,ijkEm(1)))**2 + &
                     (xyzRef(ZCOORD,ijkE(1))-xyzRef(ZCOORD,ijkEm(1)))**2)
        ds(2) = ds(2) + &
                SQRT((xyzRef(XCOORD,ijkE(2))-xyzRef(XCOORD,ijkEm(2)))**2 + &
                     (xyzRef(YCOORD,ijkE(2))-xyzRef(YCOORD,ijkEm(2)))**2 + &
                     (xyzRef(ZCOORD,ijkE(2))-xyzRef(ZCOORD,ijkEm(2)))**2)
        sum12 = .false.
      ENDIF
      ds(3) = ds(3) + &
              SQRT((xyzRef(XCOORD,ijkE(3))-xyzRef(XCOORD,ijkEm(3)))**2 + &
                   (xyzRef(YCOORD,ijkE(3))-xyzRef(YCOORD,ijkEm(3)))**2 + &
                   (xyzRef(ZCOORD,ijkE(3))-xyzRef(ZCOORD,ijkEm(3)))**2)
      ds(4) = ds(4) + &
              SQRT((xyzRef(XCOORD,ijkE(4))-xyzRef(XCOORD,ijkEm(4)))**2 + &
                   (xyzRef(YCOORD,ijkE(4))-xyzRef(YCOORD,ijkEm(4)))**2 + &
                   (xyzRef(ZCOORD,ijkE(4))-xyzRef(ZCOORD,ijkEm(4)))**2)
      s(:)  = ds(:)/arcLen(:)
      e1(:) = xyz(:,ijkE(1)) - xyzRef(:,ijkE(1))
      e2(:) = xyz(:,ijkE(2)) - xyzRef(:,ijkE(2))
      e3(:) = xyz(:,ijkE(3)) - xyzRef(:,ijkE(3))
      e4(:) = xyz(:,ijkE(4)) - xyzRef(:,ijkE(4))
      CALL RFLO_Tfint2d( s(1),s(2),s(3),s(4),e1,e2,e3,e4,p1,p2,p3,p4,dN )
      xyz(:,ijkN) = dN(:) + xyzRef(:,ijkN)
    ENDDO  ! l1
  ENDDO    ! l2

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_MoveGridCurvedPatch


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModMoveGridUtil

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModMoveGridUtil.F90,v $
! Revision 1.5  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/24 00:54:35  wasistho
! fixed loop indexing in patch%arclen1,2
!
! Revision 1.2  2006/03/18 08:19:04  wasistho
! added moveGridCurvedPatch
!
! Revision 1.1  2006/03/17 06:39:38  wasistho
! initial import
!
!
!
! ******************************************************************************








