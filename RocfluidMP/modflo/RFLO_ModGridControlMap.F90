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
! Purpose: Suite for routines computing gradients of grid control map and 
!          grid control functions.
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModGridControlMap.F90,v 1.9 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModGridControlMap

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
  PUBLIC :: RFLO_GridControlGrad3D, &
            RFLO_GridControlFunc3D, &
            RFLO_GridPhysGrad3D, &
            RFLO_GridControlGrad2D, &
            RFLO_GridControlFunc2D, &
            RFLO_GridPhysGrad2D, &
            RFLO_GridControlMap3D, &
            RFLO_GridControlMap2D

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModGridControlMap.F90,v $ $Revision: 1.9 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: compute gradient of grid control map w.r.t computational space grid.
!
! Description: the method used is second order finite difference.
!
! Input: region = grid data of current region
!
! Output: stui,stuj,stuk    = 1st order derivatives 
!         stuii,stujj,stukk = 2nd order derivatives
!         stuij,stuik,stujk = mixed derivatives 
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlGrad3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
                            
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompI, RFLO_FinDiffCompJ, &
                                       RFLO_FinDiffCompK

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: iLev,idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,iNOff,ijNOff,ndum
  REAL(RFREAL), POINTER :: stu(:,:), stui(:,:), stuj(:,:), stuk(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridControlGrad3D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions, allocate temporary storage ----------------------------------

  iLev = 1
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ndum = region%nDumCells

! get pointers ----------------------------------------------------------------

  stu   => region%levels(1)%grid%stu
  stui  => region%levels(1)%grid%stui
  stuj  => region%levels(1)%grid%stuj
  stuk  => region%levels(1)%grid%stuk

! stui ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompI( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stu,stui )

! stuj ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJ( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stu,stuj )

! stuk ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompK( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stu,stuk )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridControlGrad3D


!******************************************************************************
!
! Purpose: compute gradient of physical space grid w.r.t computational space 
!          grid.
!
! Description: the method used is second order finite difference.
!
! Input: region = grid data of current region
!
! Output: stui,stuj,stuk    = 1st order derivatives 
!         stuii,stujj,stukk = 2nd order derivatives
!         stuij,stuik,stujk = mixed derivatives 
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridPhysGrad3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
                            
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompI, RFLO_FinDiffCompJ, &
                                       RFLO_FinDiffCompK

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: iLev,idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,iNOff,ijNOff,ndum
  REAL(RFREAL), POINTER :: xyz(:,:), stui(:,:), stuj(:,:), stuk(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridPhysGrad3D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions, allocate temporary storage ----------------------------------

  iLev = 1
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ndum = region%nDumCells

! get pointers ----------------------------------------------------------------

  xyz   => region%levels(1)%gridOld%xyz
  stui  => region%levels(1)%grid%stui
  stuj  => region%levels(1)%grid%stuj
  stuk  => region%levels(1)%grid%stuk

  stui = 0._RFREAL
  stuj = 0._RFREAL
  stuk = 0._RFREAL

! stui ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompI( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,xyz,stui )

! stuj ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJ( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,xyz,stuj )

! stuk ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompK( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,xyz,stuk )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridPhysGrad3D


!******************************************************************************
!
! Purpose: compute grid control function.
!
! Description: Equation 4.65 in Handbook of Grid Generation, by
!              Joe F. Thompson, Bharat K. Soni, and Nigel P. Weatherill,
!              CRC Press 1999, ISBN 0-8493-2687-7.
!
! Input: region = grid data of current region
!
! Output: region%levels(1)%grid%pmat = Pmatrix of grid control function
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlFunc3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
                            
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompI, RFLO_FinDiffCompJ, &
                                       RFLO_FinDiffCompK, RFLO_FinDiffCompII, &
                                       RFLO_FinDiffCompJJ, RFLO_FinDiffCompKK

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ijkN, XCO, YCO, ZCO, ibn, ien, errFl
  INTEGER :: idnbeg,idnend, jdnbeg,jdnend, kdnbeg,kdnend, iNOff,ijNOff, ndum
  REAL(RFREAL) :: determ, rndet, a11, a12, a13, a21, a22, a23, a31, a32, a33
  REAL(RFREAL), POINTER :: stu(:,:), stui(:,:), stuj(:,:), stuk(:,:)
  REAL(RFREAL), POINTER :: stuii(:,:), stujj(:,:), stukk(:,:)
  REAL(RFREAL), POINTER :: stuij(:,:), stuik(:,:), stujk(:,:), pmat(:,:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_grid) , POINTER  :: grid

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridControlFunc3D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions, allocate temporary storage ----------------------------------

  XCO  = XCOORD
  YCO  = YCOORD
  ZCO  = ZCOORD

  iLev = 1
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  ndum = region%nDumCells

! allocate memory -------------------------------------------------------------

  grid => region%levels(iLev)%grid

  ALLOCATE( grid%stuii(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( grid%stujj(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( grid%stukk(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( grid%stuij(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( grid%stuik(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
  ALLOCATE( grid%stujk(    3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88

! get pointers ----------------------------------------------------------------

  stu   => region%levels(iLev)%grid%stu
  stui  => region%levels(iLev)%grid%stui
  stuj  => region%levels(iLev)%grid%stuj
  stuk  => region%levels(iLev)%grid%stuk
  stuii => region%levels(iLev)%grid%stuii
  stujj => region%levels(iLev)%grid%stujj
  stukk => region%levels(iLev)%grid%stukk
  stuij => region%levels(iLev)%grid%stuij
  stuik => region%levels(iLev)%grid%stuik
  stujk => region%levels(iLev)%grid%stujk
  pmat  => region%levels(iLev)%grid%pmat

! stuii -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompII( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                           iNOff,ijNOff,XCOORD,ZCOORD,stu,stuii )

! stujj -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompJJ( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                           iNOff,ijNOff,XCOORD,ZCOORD,stu,stujj )

! stukk -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompKK( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                           iNOff,ijNOff,XCOORD,ZCOORD,stu,stukk )

! stuij -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompJ( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stui,stuij )

! stuik -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompK( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stui,stuik )

! stujk -----------------------------------------------------------------------

  CALL RFLO_FinDiffCompK( idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,ndum, &
                          iNOff,ijNOff,XCOORD,ZCOORD,stuj,stujk )

! compute control functions ---------------------------------------------------

  DO k=kdnbeg,kdnend
    DO j=jdnbeg,jdnend
      DO i=idnbeg,idnend
        ijkN   = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        determ = stui(XCO,ijkN)*(stuj(YCO,ijkN)*stuk(ZCO,ijkN)- &
                                 stuk(YCO,ijkN)*stuj(ZCO,ijkN))+ &
                 stuj(XCO,ijkN)*(stuk(YCO,ijkN)*stui(ZCO,ijkN)- &
                                 stui(YCO,ijkN)*stuk(ZCO,ijkN))+ &
                 stuk(XCO,ijkN)*(stui(YCO,ijkN)*stuj(ZCO,ijkN)- &
                                 stuj(YCO,ijkN)*stui(ZCO,ijkN))
        rndet  = -1._RFREAL/determ
        a11    = stuj(YCO,ijkN)*stuk(ZCO,ijkN)-stuk(YCO,ijkN)*stuj(ZCO,ijkN)
        a12    = stuk(XCO,ijkN)*stuj(ZCO,ijkN)-stuj(XCO,ijkN)*stuk(ZCO,ijkN)
        a13    = stuj(XCO,ijkN)*stuk(YCO,ijkN)-stuk(XCO,ijkN)*stuj(YCO,ijkN)
        a21    = stuk(YCO,ijkN)*stui(ZCO,ijkN)-stui(YCO,ijkN)*stuk(ZCO,ijkN)
        a22    = stui(XCO,ijkN)*stuk(ZCO,ijkN)-stuk(XCO,ijkN)*stui(ZCO,ijkN)
        a23    = stuk(XCO,ijkN)*stui(YCO,ijkN)-stui(XCO,ijkN)*stuk(YCO,ijkN)
        a31    = stui(YCO,ijkN)*stuj(ZCO,ijkN)-stuj(YCO,ijkN)*stui(ZCO,ijkN)
        a32    = stuj(XCO,ijkN)*stui(ZCO,ijkN)-stui(XCO,ijkN)*stuj(ZCO,ijkN)
        a33    = stui(XCO,ijkN)*stuj(YCO,ijkN)-stuj(XCO,ijkN)*stui(YCO,ijkN)

        pmat(XCO,1,ijkN) = (stuii(XCO,ijkN)*a11+ &
                            stuii(YCO,ijkN)*a12+ &
                            stuii(ZCO,ijkN)*a13)*rndet
        pmat(YCO,1,ijkN) = (stuii(XCO,ijkN)*a21+ &
                            stuii(YCO,ijkN)*a22+ &
                            stuii(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,1,ijkN) = (stuii(XCO,ijkN)*a31+ &
                            stuii(YCO,ijkN)*a32+ &
                            stuii(ZCO,ijkN)*a33)*rndet

        pmat(XCO,2,ijkN) = (stuij(XCO,ijkN)*a11+ &
                            stuij(YCO,ijkN)*a12+ &
                            stuij(ZCO,ijkN)*a13)*rndet
        pmat(YCO,2,ijkN) = (stuij(XCO,ijkN)*a21+ &
                            stuij(YCO,ijkN)*a22+ &
                            stuij(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,2,ijkN) = (stuij(XCO,ijkN)*a31+ &
                            stuij(YCO,ijkN)*a32+ &
                            stuij(ZCO,ijkN)*a33)*rndet

        pmat(XCO,3,ijkN) = (stuik(XCO,ijkN)*a11+ &
                            stuik(YCO,ijkN)*a12+ &
                            stuik(ZCO,ijkN)*a13)*rndet
        pmat(YCO,3,ijkN) = (stuik(XCO,ijkN)*a21+ &
                            stuik(YCO,ijkN)*a22+ &
                            stuik(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,3,ijkN) = (stuik(XCO,ijkN)*a31+ &
                            stuik(YCO,ijkN)*a32+ &
                            stuik(ZCO,ijkN)*a33)*rndet

        pmat(XCO,4,ijkN) = (stujj(XCO,ijkN)*a11+ &
                            stujj(YCO,ijkN)*a12+ &
                            stujj(ZCO,ijkN)*a13)*rndet
        pmat(YCO,4,ijkN) = (stujj(XCO,ijkN)*a21+ &
                            stujj(YCO,ijkN)*a22+ &
                            stujj(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,4,ijkN) = (stujj(XCO,ijkN)*a31+ &
                            stujj(YCO,ijkN)*a32+ &
                            stujj(ZCO,ijkN)*a33)*rndet

        pmat(XCO,5,ijkN) = (stujk(XCO,ijkN)*a11+ &
                            stujk(YCO,ijkN)*a12+ &
                            stujk(ZCO,ijkN)*a13)*rndet
        pmat(YCO,5,ijkN) = (stujk(XCO,ijkN)*a21+ &
                            stujk(YCO,ijkN)*a22+ &
                            stujk(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,5,ijkN) = (stujk(XCO,ijkN)*a31+ &
                            stujk(YCO,ijkN)*a32+ &
                            stujk(ZCO,ijkN)*a33)*rndet

        pmat(XCO,6,ijkN) = (stukk(XCO,ijkN)*a11+ &
                            stukk(YCO,ijkN)*a12+ &
                            stukk(ZCO,ijkN)*a13)*rndet
        pmat(YCO,6,ijkN) = (stukk(XCO,ijkN)*a21+ &
                            stukk(YCO,ijkN)*a22+ &
                            stukk(ZCO,ijkN)*a23)*rndet
        pmat(ZCO,6,ijkN) = (stukk(XCO,ijkN)*a31+ &
                            stukk(YCO,ijkN)*a32+ &
                            stukk(ZCO,ijkN)*a33)*rndet
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! deallocate memory -----------------------------------------------------------

  DEALLOCATE( grid%stuii, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( grid%stujj, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( grid%stukk, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( grid%stuij, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( grid%stuik, stat=errFl ); IF (errFl>0) GOTO 99
  DEALLOCATE( grid%stujk, stat=errFl ); IF (errFl>0) GOTO 99

  GOTO 999

! finalize --------------------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

99   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

999  CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridControlFunc3D


!******************************************************************************
!
! Purpose: compute gradient of grid control map w.r.t computational space grid
!          for patches
!
! Description: the method used is second order finite difference.
!
! Input: region = grid data of current region
!
! Output: sti,stj   = 1st order derivatives 
!         stii,stjj = 2nd order derivatives
!         stij      = mixed derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlGrad2D( region,patch,iPatch )
                            
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompIs, RFLO_FinDiffCompJs

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... local variables
  INTEGER :: h1, h2
  REAL(RFREAL), POINTER :: st(:,:,:), sti(:,:,:), stj(:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridControlGrad2D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions --------------------------------------------------------------

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  st   => patch%st
  sti  => patch%sti
  stj  => patch%stj

! sti -------------------------------------------------------------------------

  CALL RFLO_FinDiffCompIs( h1,h2,1,2,st,sti )

! stj -------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJs( h1,h2,1,2,st,stj )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridControlGrad2D


!******************************************************************************
!
! Purpose: compute gradient of physical space grid w.r.t computational space 
!          grid for patches.
!
! Description: the method used is second order finite difference.
!
! Input: region = grid data of current region
!
! Output: sti,stj   = 1st order derivatives 
!         stii,stjj = 2nd order derivatives
!         stij      = mixed derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridPhysGrad2D( region,patch )
                            
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompIs, RFLO_FinDiffCompJs

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, h1, h2, lbound, iNOff, ijNOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ijkN, ng1, ng2
  REAL(RFREAL), POINTER :: st(:,:,:), sti(:,:,:), stj(:,:,:), xyz(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridPhysGrad2D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions --------------------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  xyz  => region%levels(iLev)%grid%xyz
  st   => patch%st
  sti  => patch%sti
  stj  => patch%stj

! copy xyz to st --------------------------------------------------------------

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkN = IndIJK(i,j,k,iNOff,ijNOff)
        IF      (lbound==1 .OR. lbound==2) THEN
          ng1 = j - jbeg + 1
          ng2 = k - kbeg + 1
        ELSE IF (lbound==3 .OR. lbound==4) THEN
          ng1 = k - kbeg + 1
          ng2 = i - ibeg + 1
        ELSE IF (lbound==5 .OR. lbound==6) THEN
          ng1 = i - ibeg + 1
          ng2 = j - jbeg + 1
        ENDIF
        st(1,ng1,ng2) = xyz(XCOORD,ijkN)
        st(2,ng1,ng2) = xyz(YCOORD,ijkN)
        st(3,ng1,ng2) = xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

! sti -------------------------------------------------------------------------

  CALL RFLO_FinDiffCompIs( h1,h2,1,3,st,sti )

! stj -------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJs( h1,h2,1,3,st,stj )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridPhysGrad2D


!******************************************************************************
!
! Purpose: compute grid control function.
!
! Description: Equation 4.13 in Handbook of Grid Generation, by
!              Joe F. Thompson, Bharat K. Soni, and Nigel P. Weatherill,
!              CRC Press 1999, ISBN 0-8493-2687-7.
!
! Input: region = grid data of current region
!
! Output: patch%pfun = Pmatrix of grid control function
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlFunc2D( region,patch,iPatch )
                            
  USE RFLO_ModFiniteDifference, ONLY : RFLO_FinDiffCompIs, &
                                       RFLO_FinDiffCompJs, &
                                       RFLO_FinDiffCompIIs, &
                                       RFLO_FinDiffCompJJs

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: iLev, lbound, h1, h2
  REAL(RFREAL) :: determ, rndet, a11, a12, a21, a22
  REAL(RFREAL), POINTER :: st(:,:,:), sti(:,:,:), stj(:,:,:)
  REAL(RFREAL), POINTER :: stii(:,:,:), stjj(:,:,:), stij(:,:,:), pfun(:,:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridControlFunc2D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions --------------------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  st   => patch%st
  sti  => patch%sti
  stj  => patch%stj
  stii => patch%stii
  stjj => patch%stjj
  stij => patch%stij
  pfun => patch%pfun

! stii ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompIIs( h1,h2,1,3,st,stii )

! stjj ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJJs( h1,h2,1,3,st,stjj )

! stij ------------------------------------------------------------------------

  CALL RFLO_FinDiffCompJs( h1,h2,1,3,sti,stij )

! compute control functions ---------------------------------------------------

  DO j=1,h2
    DO i=1,h1
      determ = sti(1,i,j)*stj(2,i,j)-stj(1,i,j)*sti(2,i,j)
      rndet  = -1._RFREAL/determ
      a11    =  stj(2,i,j)
      a12    = -stj(1,i,j)
      a21    = -sti(2,i,j)
      a22    =  sti(1,i,j)

      pfun(1,1,i,j) = (stii(1,i,j)*a11 + stii(2,i,j)*a12)*rndet
      pfun(2,1,i,j) = (stii(1,i,j)*a12 + stii(2,i,j)*a22)*rndet

      pfun(1,2,i,j) = (stij(1,i,j)*a11 + stij(2,i,j)*a12)*rndet
      pfun(2,2,i,j) = (stij(1,i,j)*a12 + stij(2,i,j)*a22)*rndet

      pfun(1,3,i,j) = (stjj(1,i,j)*a11 + stjj(2,i,j)*a12)*rndet
      pfun(2,3,i,j) = (stjj(1,i,j)*a12 + stjj(2,i,j)*a22)*rndet
    ENDDO  ! i
  ENDDO    ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridControlFunc2D


!******************************************************************************
!
! Purpose: redistribute the control parameter grid based on the normalized arc-
!          length of block edges.
!
! Description: the method used is bi-linear interpolation.
!
! Input: region = grid data of current region
!
! Output: stu = new parameter grid.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlMap3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                            RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset

  USE RFLO_ModExtrapolation, ONLY : RFLO_ExtrapRegDummyNode 

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iter

! ... local variables
  INTEGER :: iLev, ibeg, iend, jbeg, jend, kbeg, kend, ndum
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijkN, imjkN, ijmkN, ijkmN, iNOff, ijNOff, errorFlag
  INTEGER :: XCO, YCO, ZCO

  REAL(RFREAL) :: phii, phii1, phij, phij1, phik, phik1, dsi, dx,dy,dz, resid
  REAL(RFREAL), POINTER :: xyzIni(:,:), stu(:,:), stuOld(:,:)
  REAL(RFREAL), ALLOCATABLE :: dsj(:), dsk(:,:)
  REAL(RFREAL), ALLOCATABLE :: arclen12(:,:), arclen34(:,:), arclen56(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GridControlMap3D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions, allocate temporary storage ----------------------------------

  XCO = XCOORD
  YCO = YCOORD
  ZCO = ZCOORD

  ndum = region%nDumCells

  iLev = 1
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ALLOCATE( dsj(ipnbeg:ipnend)              ,stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( dsk(ipnbeg:ipnend,jpnbeg:jpnend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( arclen12(jpnbeg:jpnend,kpnbeg:kpnend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( arclen34(kpnbeg:kpnend,ipnbeg:ipnend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

  ALLOCATE( arclen56(ipnbeg:ipnend,jpnbeg:jpnend),stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

! get pointers ----------------------------------------------------------------

  xyzIni => region%levels(1)%gridOld%xyzOld
  stu    => region%levels(1)%grid%stu
  stuOld => region%levels(1)%grid%stuOld

! initialize stuOld -----------------------------------------------------------

  stuOld   = 0._RFREAL
  arclen12 = 0._RFREAL
  arclen34 = 0._RFREAL
  arclen56 = 0._RFREAL

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg+1,ipnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)

        arclen12(j,k) = arclen12(j,k) + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,imjkN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,imjkN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,imjkN))**2)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO i=ipnbeg,ipnend
    DO k=kpnbeg,kpnend
      DO j=jpnbeg+1,jpnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)

        arclen34(k,i) = arclen34(k,i) + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijmkN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijmkN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijmkN))**2)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO j=jpnbeg,jpnend
    DO i=ipnbeg,ipnend
      DO k=kpnbeg+1,kpnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        arclen56(i,j) = arclen56(i,j) + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkmN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkmN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkmN))**2)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      dsi  = 0._RFREAL
      DO i=ipnbeg+1,ipnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)

        dsi      = dsi      + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,imjkN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,imjkN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,imjkN))**2)

        stuOld(XCO,ijkN) = dsi/arcLen12(j,k)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO k=kpnbeg,kpnend
    dsj(:) = 0._RFREAL
    DO j=jpnbeg+1,jpnend
      DO i=ipnbeg,ipnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)

        dsj(i)   = dsj(i)   + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijmkN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijmkN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijmkN))**2)

        stuOld(YCO,ijkN) = dsj(i)/arcLen34(k,i)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  dsk(:,:) = 0._RFREAL
  DO k=kpnbeg+1,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

        dsk(i,j) = dsk(i,j) + &
                     SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkmN))**2 + &
                          (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkmN))**2 + &
                          (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkmN))**2)

        stuOld(ZCO,ijkN) = dsk(i,j)/arcLen56(i,j)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! iterations on control-parameter coordinates --------------------------------

  stu = stuOld

  DO iter=1,5

! - interpolate parameter coordinates inside region --------------------------

    dsk(:,:) = 0._RFREAL
    DO k=kpnbeg+1,kpnend-1
      dsj(:) = 0._RFREAL
      DO j=jpnbeg+1,jpnend-1
        dsi  = 0._RFREAL
        DO i=ipnbeg+1,ipnend-1
          ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
          ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
          ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)

          dsi      = dsi      + stuOld(XCOORD,ijkN)-stuOld(XCOORD,imjkN)
          dsj(i)   = dsj(i)   + stuOld(YCOORD,ijkN)-stuOld(YCOORD,ijmkN)
          dsk(i,j) = dsk(i,j) + stuOld(ZCOORD,ijkN)-stuOld(ZCOORD,ijkmN)
                            

          phii  = dsi
          phii1 = 1._RFREAL - phii
          phij  = dsj(i)
          phij1 = 1._RFREAL - phij
          phik  = dsk(i,j)
          phik1 = 1._RFREAL - phik

          stu(XCO,ijkN) = phij1*phik1* &
                          stuOld(XCO,IndIJK(i,jpnbeg,kpnbeg,iNOff,ijNOff)) + &
                          phij*phik* &
                          stuOld(XCO,IndIJK(i,jpnend,kpnend,iNOff,ijNOff)) + &
                          phij*phik1* &
                          stuOld(XCO,IndIJK(i,jpnend,kpnbeg,iNOff,ijNOff)) + &
                          phij1*phik* &
                          stuOld(XCO,IndIJK(i,jpnbeg,kpnend,iNOff,ijNOff))

          stu(YCO,ijkN) = phii1*phik1* &
                          stuOld(YCO,IndIJK(ipnbeg,j,kpnbeg,iNOff,ijNOff)) + &
                          phii*phik* &
                          stuOld(YCO,IndIJK(ipnend,j,kpnend,iNOff,ijNOff)) + &
                          phii*phik1* &
                          stuOld(YCO,IndIJK(ipnend,j,kpnbeg,iNOff,ijNOff)) + &
                          phii1*phik* &
                          stuOld(YCO,IndIJK(ipnbeg,j,kpnend,iNOff,ijNOff))

          stu(ZCO,ijkN) = phii1*phij1* &
                          stuOld(ZCO,IndIJK(ipnbeg,jpnbeg,k,iNOff,ijNOff)) + &
                          phii*phij* &
                          stuOld(ZCO,IndIJK(ipnend,jpnend,k,iNOff,ijNOff)) + &
                          phii*phij1* &
                          stuOld(ZCO,IndIJK(ipnend,jpnbeg,k,iNOff,ijNOff)) + &
                          phii1*phij* &
                          stuOld(ZCO,IndIJK(ipnbeg,jpnend,k,iNOff,ijNOff))
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! - residual ------------------------------------------------------------------

    resid = 0._RFREAL
    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          ijkN  = IndIJK(i,j,k,iNOff,ijNOff)
          dx    = stu(XCOORD,ijkN) - stuOld(XCOORD,ijkN)
          dy    = stu(YCOORD,ijkN) - stuOld(YCOORD,ijkN)
          dz    = stu(ZCOORD,ijkN) - stuOld(ZCOORD,ijkN)
          resid = resid + dx*dx + dy*dy +dz*dz
        ENDDO
      ENDDO
    ENDDO

    IF (global%myProcid == (global%nProcAlloc-1)/2 .AND. &
        global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,*) SOLVER_NAME//' CtrParam_Vol:region,iter,residual', &
                                     region%iRegionGlobal,iter,resid
    ENDIF

! - assign stuOld -------------------------------------------------------------

    stuOld = stu

  ENDDO  ! iter

! extrapolate stu to dummy ----------------------------------------------------

  ibeg = idnbeg
  iend = idnend
  jbeg = jdnbeg
  jend = jdnend
  kbeg = kdnbeg
  kend = kdnend

  CALL RFLO_ExtrapRegDummyNode( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                                iNOff,ijNOff,XCOORD,ZCOORD,stu )

! finalize --------------------------------------------------------------------

  DEALLOCATE( dsj,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( dsk,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( arclen12,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( arclen34,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  DEALLOCATE( arclen56,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridControlMap3D


!******************************************************************************
!
! Purpose: redistribute the control parameter grid based on the normalized arc-
!          length of patch edges.
!
! Description: the method used is Hermite cubic interpolation.
!
! Input: region = grid data of current region
!        patch  = grid data of current patch
!
! Output: st = new parameter grid.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridControlMap2D( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensPhysNodes

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: iPatch

! ... loop variables
  INTEGER :: i, j, k, iter

! ... local variables
  INTEGER :: iLev, lbound, h1, h2, ir, jr, kr, ijkN, ijkNm, m1, m2, errorFlag
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff

  REAL(RFREAL) :: arclen(4), ri, rj, dsi, phii, phii1, phij, phij1
  REAL(RFREAL) :: dx, dy, resid
  REAL(RFREAL), ALLOCATABLE :: dsj(:)
  REAL(RFREAL), POINTER :: xyzIni(:,:), st(:,:,:), stOld(:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_GridControlMap2D',&
       'RFLO_ModGridControlMap.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyzIni => region%levels(iLev)%gridOld%xyzOld
  st     => patch%st
  stOld  => patch%stOld

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

  ALLOCATE( dsj(h1), stat=errorFlag )
  global%error = errorFlag
  IF (global%error/=0) CALL ErrorStop( global,ERR_ALLOCATE,&
       __LINE__ )

! initialize stOld ------------------------------------------------------------

  stOld     = 0._RFREAL
  arclen(:) = 0._RFREAL

  IF (lbound==1 .OR. lbound==2) THEN
    IF (lbound==1) THEN
      ir = ipnbeg
    ELSE
      ir = ipnend
    ENDIF
    DO k=kpnbeg+1,kpnend
      ijkN    = IndIJK(ir ,jpnbeg,k  ,iNOff,ijNOff)
      ijkNm   = IndIJK(ir ,jpnbeg,k-1,iNOff,ijNOff)
      m1      = 1
      m2      = k-kpnbeg+1
      arclen(1) = arclen(1) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(1)
    ENDDO       ! k
    stOld(1,m1,:) = 0._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(1)

    DO k=kpnbeg+1,kpnend
      ijkN    = IndIJK(ir ,jpnend,k  ,iNOff,ijNOff)
      ijkNm   = IndIJK(ir ,jpnend,k-1,iNOff,ijNOff)
      m1      = jpnend-jpnbeg+1
      m2      = k-kpnbeg+1
      arclen(2) = arclen(2) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(2)
    ENDDO       ! k
    stOld(1,m1,:) = 1._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(2)

    DO j=jpnbeg+1,jpnend
      ijkN    = IndIJK(ir ,j  ,kpnbeg,iNOff,ijNOff)
      ijkNm   = IndIJK(ir ,j-1,kpnbeg,iNOff,ijNOff)
      m1      = j-jpnbeg+1
      m2      = 1
      arclen(3) = arclen(3) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(3)
    ENDDO       ! j
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(3)
    stOld(2,:,m2) = 0._RFREAL

    DO j=jpnbeg+1,jpnend
      ijkN    = IndIJK(ir ,j  ,kpnbeg,iNOff,ijNOff)
      ijkNm   = IndIJK(ir ,j-1,kpnbeg,iNOff,ijNOff)
      m1      = j-jpnbeg+1
      m2      = kpnend-kpnbeg+1
      arclen(4) = arclen(4) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(4)
    ENDDO       ! j
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(4)
    stOld(2,:,m2) = 1._RFREAL

  ELSEIF (lbound==3 .OR. lbound==4) THEN
    IF (lbound==3) THEN
      jr = jpnbeg
    ELSE
      jr = jpnend
    ENDIF
    DO i=ipnbeg+1,ipnend
      ijkN    = IndIJK(i  ,jr ,kpnbeg,iNOff,ijNOff)
      ijkNm   = IndIJK(i-1,jr ,kpnbeg,iNOff,ijNOff)
      m1      = 1
      m2      = i-ipnbeg+1
      arclen(1) = arclen(1) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(1)
    ENDDO       ! i
    stOld(1,m1,:) = 0._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(1)

    DO i=ipnbeg+1,ipnend
      ijkN    = IndIJK(i  ,jr ,kpnend,iNOff,ijNOff)
      ijkNm   = IndIJK(i-1,jr ,kpnend,iNOff,ijNOff)
      m1      = kpnend-kpnbeg+1
      m2      = i-ipnbeg+1
      arclen(2) = arclen(2) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(2)
    ENDDO       ! i
    stOld(1,m1,:) = 1._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(2)

    DO k=kpnbeg+1,kpnend
      ijkN    = IndIJK(ipnbeg,jr ,k  ,iNOff,ijNOff)
      ijkNm   = IndIJK(ipnbeg,jr ,k-1,iNOff,ijNOff)
      m1      = k-kpnbeg+1
      m2      = 1
      arclen(3) = arclen(3) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(3)
    ENDDO       ! k
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(3)
    stOld(2,:,m2) = 0._RFREAL

    DO k=kpnbeg+1,kpnend
      ijkN    = IndIJK(ipnend,jr ,k  ,iNOff,ijNOff)
      ijkNm   = IndIJK(ipnend,jr ,k-1,iNOff,ijNOff)
      m1      = k-kpnbeg+1
      m2      = ipnend-ipnbeg+1
      arclen(4) = arclen(4) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(4)
    ENDDO       ! k
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(4)
    stOld(2,:,m2) = 1._RFREAL

  ELSEIF (lbound==5 .OR. lbound==6) THEN
    IF (lbound==5) THEN
      kr = kpnbeg
    ELSE
      kr = kpnend
    ENDIF
    DO j=jpnbeg+1,jpnend
      ijkN    = IndIJK(ipnbeg,j  ,kr ,iNOff,ijNOff)
      ijkNm   = IndIJK(ipnbeg,j-1,kr ,iNOff,ijNOff)
      m1      = 1
      m2      = j-jpnbeg+1
      arclen(1) = arclen(1) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(1)
    ENDDO       ! j
    stOld(1,m1,:) = 0._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(1)

    DO j=jpnbeg+1,jpnend
      ijkN    = IndIJK(ipnend,j  ,kr ,iNOff,ijNOff)
      ijkNm   = IndIJK(ipnend,j-1,kr ,iNOff,ijNOff)
      m1      = ipnend-ipnbeg+1
      m2      = j-jpnbeg+1
      arclen(2) = arclen(2) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(2,m1,m2) = arclen(2)
    ENDDO       ! j
    stOld(1,m1,:) = 1._RFREAL
    stOld(2,m1,:) = stOld(2,m1,:)/arclen(2)

    DO i=ipnbeg+1,ipnend
      ijkN    = IndIJK(i  ,jpnbeg,kr ,iNOff,ijNOff)
      ijkNm   = IndIJK(i-1,jpnbeg,kr ,iNOff,ijNOff)
      m1      = i-ipnbeg+1
      m2      = 1
      arclen(3) = arclen(3) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(3)
    ENDDO       ! i
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(3)
    stOld(2,:,m2) = 0._RFREAL

    DO i=ipnbeg+1,ipnend
      ijkN    = IndIJK(i  ,jpnend,kr ,iNOff,ijNOff)
      ijkNm   = IndIJK(i-1,jpnend,kr ,iNOff,ijNOff)
      m1      = i-ipnbeg+1
      m2      = jpnend-jpnbeg+1
      arclen(4) = arclen(4) + &
                  SQRT((xyzIni(XCOORD,ijkN)-xyzIni(XCOORD,ijkNm))**2 + &
                       (xyzIni(YCOORD,ijkN)-xyzIni(YCOORD,ijkNm))**2 + &
                       (xyzIni(ZCOORD,ijkN)-xyzIni(ZCOORD,ijkNm))**2)
      stOld(1,m1,m2) = arclen(4)
    ENDDO       ! i
    stOld(1,:,m2) = stOld(1,:,m2)/arclen(4)
    stOld(2,:,m2) = 1._RFREAL

  ENDIF         ! lbound

  DO j=1,h2
    DO i=2,h1-1
      ri = REAL(i-1)/REAL(h1-1)
      stOld(2,i,j) = ri*stOld(2,h1,j) + (1._RFREAL-ri)*stOld(2,1,j)
    ENDDO  ! i
  ENDDO    ! j

  DO j=2,h2-1
    DO i=1,h1
      rj = REAL(j-1)/REAL(h2-1)
      stOld(1,i,j) = rj*stOld(1,i,h2) + (1._RFREAL-rj)*stOld(1,i,1)
    ENDDO  ! i
  ENDDO    ! j

! iterations on control-parameter coordinates --------------------------------

  DO iter=1,5

! - interpolate parameter coordinates inside region --------------------------

    st = stOld

    DO j=1,h2
      dsi  = 0._RFREAL
      DO i=2,h1
        dsi   = dsi + stOld(1,i,j)-stOld(1,i-1,j)

        phii  = (3._RFREAL-2._RFREAL*dsi)*dsi*dsi
        phii1 = 1._RFREAL - phii
        st(2,i,j) = phii*stOld(2,h1,j) + phii1*stOld(2,1,j)
      ENDDO   ! i
    ENDDO     ! j

    dsj(:) = 0._RFREAL
    DO j=2,h2
      DO i=1,h1
        dsj(i)   = dsj(i) + stOld(2,i,j)-stOld(2,i,j-1)

        phij  = (3._RFREAL-2._RFREAL*dsj(i))*dsj(i)*dsj(i)
        phij1 = 1._RFREAL - phij

        st(1,i,j) = phij*stOld(1,i,h2) + phij1*stOld(1,i,1)
      ENDDO   ! i
    ENDDO     ! j

! - residual ------------------------------------------------------------------

    resid = 0._RFREAL
    DO j=1,h2
      DO i=1,h1
        dx    = st(1,i,j) - stOld(1,i,j)
        dy    = st(2,i,j) - stOld(2,i,j)
        resid = resid + dx*dx + dy*dy
      ENDDO
    ENDDO

    IF (global%myProcid == (global%nProcAlloc-1)/2 .AND. &
        global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,*) SOLVER_NAME//' CtrParam_Surf:patch,iter,residual', &
                                     region%iRegionGlobal,iPatch,iter,resid
    ENDIF

! - assign stOld --------------------------------------------------------------

    stOld = st

  ENDDO  ! iter

! finalize --------------------------------------------------------------------

  DEALLOCATE( dsj,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,&
       __LINE__ )

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GridControlMap2D


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModGridControlMap


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModGridControlMap.F90,v $
! Revision 1.9  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/03/02 00:23:59  wasistho
! prepared elliptic pde grid motion
!
! Revision 1.6  2006/02/24 21:45:26  wasistho
! bug fixed in computing stuOld in controlMap3
!
! Revision 1.5  2006/02/09 07:41:49  wasistho
! bug fixed: computation of stuj
!
! Revision 1.4  2006/02/09 00:24:13  wasistho
! bug fixed allocation of dsj in Map2
!
! Revision 1.3  2006/02/08 07:51:26  wasistho
! added iPatch in controlMap2
!
! Revision 1.2  2005/12/07 08:46:04  wasistho
! added stuff for surface mesh motion EPDE
!
! Revision 1.1  2005/12/03 09:39:47  wasistho
! initial import
!
!
!
! ******************************************************************************














