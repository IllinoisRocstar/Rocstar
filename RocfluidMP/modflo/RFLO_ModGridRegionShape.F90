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
! Purpose: Suite for region shape related routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModGridRegionShape.F90,v 1.3 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModGridRegionShape

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
  PUBLIC :: RFLO_GridFlatPatch, &
            RFLO_FindFunkyBlocks

! private : RFLO_BlockSkewedCell
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModGridRegionShape.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: identify flat patches
!
! Description: none.
!
! Input: region = data of current region, grid movements
!        patch  = current patch.
!
! Output: patch%bndFlat get value true or false
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridFlatPatch( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE RFLO_ModVectorTensor, ONLY : RFLO_NormCrossProd

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, lbound, ihlf, jhlf, khlf
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff

  REAL(RFREAL) :: vmag, eps
  REAL(RFREAL) :: v11(XCOORD:ZCOORD), v12(XCOORD:ZCOORD)
  REAL(RFREAL) :: v21(XCOORD:ZCOORD), v22(XCOORD:ZCOORD)
  REAL(RFREAL) :: v31(XCOORD:ZCOORD), v32(XCOORD:ZCOORD)
  REAL(RFREAL) :: v41(XCOORD:ZCOORD), v42(XCOORD:ZCOORD)
  REAL(RFREAL) ::  v1(XCOORD:ZCOORD),  v2(XCOORD:ZCOORD)
  REAL(RFREAL) ::  v3(XCOORD:ZCOORD),  v4(XCOORD:ZCOORD), vec(XCOORD:ZCOORD)
  REAL(RFREAL), POINTER :: xyz(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GridFlatPatch',&
       'RFLO_ModGridRegionShape.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz => region%levels(iLev)%grid%xyz

  ihlf = ibeg+(iend-ibeg)/2
  jhlf = jbeg+(jend-jbeg)/2
  khlf = kbeg+(kend-kbeg)/2

! default values --------------------------------------------------------------

  patch%bndFlat = .FALSE.
  patch%dirFlat = -1

! four pair vectors with consistent orientation -------------------------------

  IF (lbound==1 .OR. lbound==2) THEN
    v11(:) = xyz(:,IndIJK(ibeg,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v12(:) = xyz(:,IndIJK(ibeg,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v21(:) = xyz(:,IndIJK(ibeg,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff))
    v22(:) = xyz(:,IndIJK(ibeg,jhlf,kend,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff))
    v31(:) = xyz(:,IndIJK(ibeg,jhlf,kend,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kend,iNOff,ijNOff))
    v32(:) = xyz(:,IndIJK(ibeg,jend,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kend,iNOff,ijNOff))
    v41(:) = xyz(:,IndIJK(ibeg,jend,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff))
    v42(:) = xyz(:,IndIJK(ibeg,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff))
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    v11(:) = xyz(:,IndIJK(ihlf,jbeg,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v12(:) = xyz(:,IndIJK(ibeg,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v21(:) = xyz(:,IndIJK(ibeg,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff))
    v22(:) = xyz(:,IndIJK(ihlf,jbeg,kend,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kend,iNOff,ijNOff))
    v31(:) = xyz(:,IndIJK(ihlf,jbeg,kend,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kend,iNOff,ijNOff))
    v32(:) = xyz(:,IndIJK(iend,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kend,iNOff,ijNOff))
    v41(:) = xyz(:,IndIJK(iend,jbeg,khlf,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff))
    v42(:) = xyz(:,IndIJK(ihlf,jbeg,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff))
  ELSEIF (lbound==5 .OR. lbound==6) THEN
    v11(:) = xyz(:,IndIJK(ihlf,jbeg,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v12(:) = xyz(:,IndIJK(ibeg,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jbeg,kbeg,iNOff,ijNOff))
    v21(:) = xyz(:,IndIJK(ibeg,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff))
    v22(:) = xyz(:,IndIJK(ihlf,jend,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(ibeg,jend,kbeg,iNOff,ijNOff))
    v31(:) = xyz(:,IndIJK(ihlf,jend,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jend,kbeg,iNOff,ijNOff))
    v32(:) = xyz(:,IndIJK(iend,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jend,kbeg,iNOff,ijNOff))
    v41(:) = xyz(:,IndIJK(iend,jhlf,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff))
    v42(:) = xyz(:,IndIJK(ihlf,jbeg,kbeg,iNOff,ijNOff)) - &
             xyz(:,IndIJK(iend,jbeg,kbeg,iNOff,ijNOff))
  ENDIF

! unit vectors orthogonal to the patch at four patch corners ------------------

  CALL RFLO_NormCrossProd( v11,v12,v1 )
  CALL RFLO_NormCrossProd( v21,v22,v2 )
  CALL RFLO_NormCrossProd( v31,v32,v3 )
  CALL RFLO_NormCrossProd( v41,v42,v4 )

! their total length should equals 4. for flat patch --------------------------

  eps  = 1.E-3_RFREAL 
  vec  = v1+v2+v3+v4
  vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
  IF (ABS( vmag-4._RFREAL) < eps) THEN
    patch%bndFlat = .TRUE.
    patch%dirFlat = OFF
  ELSE
    region%levels(iLev)%grid%boundFlat(lbound) = .FALSE.
  ENDIF

! search for bend edges -------------------------------------------------------

  IF (lbound==1) THEN
    v1   = v12/SQRT( v12(XCOORD)**2 + v12(YCOORD)**2 + v12(ZCOORD)**2 )
    v2   = v21/SQRT( v21(XCOORD)**2 + v21(YCOORD)**2 + v21(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(1) = .FALSE.
    ENDIF

    v1   = v22/SQRT( v22(XCOORD)**2 + v22(YCOORD)**2 + v22(ZCOORD)**2 )
    v2   = v31/SQRT( v31(XCOORD)**2 + v31(YCOORD)**2 + v31(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(2) = .FALSE.
    ENDIF

    v1   = v32/SQRT( v32(XCOORD)**2 + v32(YCOORD)**2 + v32(ZCOORD)**2 )
    v2   = v41/SQRT( v41(XCOORD)**2 + v41(YCOORD)**2 + v41(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(3) = .FALSE.
    ENDIF

    v1   = v42/SQRT( v42(XCOORD)**2 + v42(YCOORD)**2 + v42(ZCOORD)**2 )
    v2   = v11/SQRT( v11(XCOORD)**2 + v11(YCOORD)**2 + v11(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(4) = .FALSE.
    ENDIF

  ELSEIF (lbound==2) THEN
    v1   = v12/SQRT( v12(XCOORD)**2 + v12(YCOORD)**2 + v12(ZCOORD)**2 )
    v2   = v21/SQRT( v21(XCOORD)**2 + v21(YCOORD)**2 + v21(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(5) = .FALSE.
    ENDIF

    v1   = v22/SQRT( v22(XCOORD)**2 + v22(YCOORD)**2 + v22(ZCOORD)**2 )
    v2   = v31/SQRT( v31(XCOORD)**2 + v31(YCOORD)**2 + v31(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(6) = .FALSE.
    ENDIF

    v1   = v32/SQRT( v32(XCOORD)**2 + v32(YCOORD)**2 + v32(ZCOORD)**2 )
    v2   = v41/SQRT( v41(XCOORD)**2 + v41(YCOORD)**2 + v41(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(7) = .FALSE.
    ENDIF

    v1   = v42/SQRT( v42(XCOORD)**2 + v42(YCOORD)**2 + v42(ZCOORD)**2 )
    v2   = v11/SQRT( v11(XCOORD)**2 + v11(YCOORD)**2 + v11(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(8) = .FALSE.
    ENDIF

  ELSEIF (lbound==3) THEN
    v1   = v11/SQRT( v11(XCOORD)**2 + v11(YCOORD)**2 + v11(ZCOORD)**2 )
    v2   = v42/SQRT( v42(XCOORD)**2 + v42(YCOORD)**2 + v42(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(9) = .FALSE.
    ENDIF

    v1   = v22/SQRT( v22(XCOORD)**2 + v22(YCOORD)**2 + v22(ZCOORD)**2 )
    v2   = v31/SQRT( v31(XCOORD)**2 + v31(YCOORD)**2 + v31(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(10) = .FALSE.
    ENDIF

  ELSEIF (lbound==4) THEN
    v1   = v11/SQRT( v11(XCOORD)**2 + v11(YCOORD)**2 + v11(ZCOORD)**2 )
    v2   = v42/SQRT( v42(XCOORD)**2 + v42(YCOORD)**2 + v42(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(11) = .FALSE. ! officially this
    ENDIF                                                 ! should be edge 12

    v1   = v22/SQRT( v22(XCOORD)**2 + v22(YCOORD)**2 + v22(ZCOORD)**2 )
    v2   = v31/SQRT( v31(XCOORD)**2 + v31(YCOORD)**2 + v31(ZCOORD)**2 )
    vec  = v1 - v2
    vmag = SQRT( vec(XCOORD)**2 + vec(YCOORD)**2 + vec(ZCOORD)**2 )
    IF (ABS( vmag-2._RFREAL) > eps) THEN
      region%levels(iLev)%grid%edgeStraight(12) = .FALSE. ! officially this
    ENDIF                                                 ! should be edge 11

  ENDIF

! search for flat direction on curved patches ---------------------------------

  IF (.NOT. patch%bndFlat) THEN
    IF (lbound==1) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(1) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(3) .EQV. .TRUE.)) THEN
        patch%dirFlat = KCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(2) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(4) .EQV. .TRUE.)) THEN 
        patch%dirFlat = JCOORD
      ENDIF
    ELSEIF (lbound==2) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(5) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(7) .EQV. .TRUE.)) THEN
        patch%dirFlat = KCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(6) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(8) .EQV. .TRUE.)) THEN 
        patch%dirFlat = JCOORD
      ENDIF
    ELSEIF (lbound==3) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(1) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(5) .EQV. .TRUE.)) THEN
        patch%dirFlat = KCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(9) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(10) .EQV. .TRUE.)) THEN 
        patch%dirFlat = ICOORD
      ENDIF
    ELSEIF (lbound==4) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(3) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(7) .EQV. .TRUE.)) THEN
        patch%dirFlat = KCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(11) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(12) .EQV. .TRUE.)) THEN 
        patch%dirFlat = ICOORD
      ENDIF
    ELSEIF (lbound==5) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(4) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(8) .EQV. .TRUE.)) THEN
        patch%dirFlat = JCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(9) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(11) .EQV. .TRUE.)) THEN 
        patch%dirFlat = ICOORD
      ENDIF
    ELSEIF (lbound==6) THEN
      IF ((region%levels(iLev)%grid%edgeStraight(2) .EQV. .TRUE.) .AND. &
          (region%levels(iLev)%grid%edgeStraight(6) .EQV. .TRUE.)) THEN
        patch%dirFlat = JCOORD
      ELSEIF ((region%levels(iLev)%grid%edgeStraight(10) .EQV. .TRUE.) .AND. &
              (region%levels(iLev)%grid%edgeStraight(12) .EQV. .TRUE.)) THEN 
        patch%dirFlat = ICOORD
      ENDIF ! edgeStraight
    ENDIF   ! lbound
  ENDIF     ! patch not flat
        
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GridFlatPatch


!******************************************************************************
!
! Purpose: identify blocks with unconventional shape, if any.
!
! Description: none.
!
! Input: regions%grid = dimensions, grid coordinates.
!
! Output: regions%blockShape = normal or funky
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FindFunkyBlocks( regions )

  USE ModInterfaces, ONLY : RFLO_CalcFaceVectors

  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'RFLO_FindFunkyBlocks',&
       'RFLO_ModGridRegionShape.F90' )

!  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
!       regions(1)%global%verbLevel > VERBOSE_NONE ) THEN
!    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Entering RFLO_FindFunkyBlocks...'
!  END IF ! global%verbLevel

! loop over all regions

  DO iReg=1,regions(1)%global%nRegions
    IF (regions(iReg)%procid==regions(iReg)%global%myProcid & ! region active and
        .AND. regions(iReg)%active==ACTIVE) THEN              ! on my processor

      CALL RFLO_CalcFaceVectors( regions(iReg) )
      CALL RFLO_BlockSkewedCell( regions(iReg) )

    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

! finalize --------------------------------------------------------------------

!  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
!       regions(1)%global%verbLevel > VERBOSE_NONE ) THEN
!    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Leaving RFLO_FindFunkyBlocks.'
!  END IF ! global%verbLevel

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE RFLO_FindFunkyBlocks


!******************************************************************************
!
! Purpose: identify if this region contains skewed cells.
!
! Description: if yes, mark it as funky block, otherwise normal block.
!
! Input: region%levels%grid = dimensions, coordinates, face vectors
!                             (current region)
!
! Output: for finest grid only.
!
!******************************************************************************

SUBROUTINE RFLO_BlockSkewedCell( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iNOff, ijNOff, ijkNode(8), face(6)
  INTEGER :: im, jm, km

  REAL(RFREAL) :: edge(3,4), dprod(8)
  REAL(RFREAL) :: dprodm, dpmax, emag, smag1, smag2
  REAL(RFREAL), POINTER :: xyz(:,:), si(:,:), sj(:,:), sk(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_BlockSkewedCell',&
       'RFLO_ModGridRegionShape.F90' )

! compute cell skewness -------------------------------------------------------

  iLev=1

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz => region%levels(iLev)%grid%xyz
  si  => region%levels(iLev)%grid%si
  sj  => region%levels(iLev)%grid%sj
  sk  => region%levels(iLev)%grid%sk

! check face vectors: cell skewness

  dprodm =  1._RFREAL
  dpmax  = -1._RFREAL

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkNode(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ijkNode(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijkNode(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkNode(4) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        ijkNode(5) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
        ijkNode(6) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
        ijkNode(7) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        ijkNode(8) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

        face(1) = ijkNode(1)
        face(2) = ijkNode(2)
        face(4) = ijkNode(3)
        face(6) = ijkNode(4)

! ----- check for skewed i-faces

        edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(2))
        edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(2))
        edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(2))

        edge(1,2) = xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(6))
        edge(2,2) = xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(6))
        edge(3,2) = xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(6))

        edge(1,3) = xyz(XCOORD,ijkNode(8))-xyz(XCOORD,ijkNode(5))
        edge(2,3) = xyz(YCOORD,ijkNode(8))-xyz(YCOORD,ijkNode(5))
        edge(3,3) = xyz(ZCOORD,ijkNode(8))-xyz(ZCOORD,ijkNode(5))

        edge(1,4) = xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(7))
        edge(2,4) = xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(7))
        edge(3,4) = xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(7))

        DO l=1,4
          emag  = SQRT( edge(1,l)*edge(1,l) + &
                        edge(2,l)*edge(2,l) + &
                        edge(3,l)*edge(3,l) )
          smag1 = SQRT( si(1,face(1))*si(1,face(1)) + &
                        si(2,face(1))*si(2,face(1)) + &
                        si(3,face(1))*si(3,face(1)) )
          smag2 = SQRT( si(1,face(2))*si(1,face(2)) + &
                        si(2,face(2))*si(2,face(2)) + &
                        si(3,face(2))*si(3,face(2)) )

          dprod(l)  =  (si(1,face(1))*edge(1,l) + &
                        si(2,face(1))*edge(2,l) + &
                        si(3,face(1))*edge(3,l))/(smag1*emag)

          dprod(l+4) = (si(1,face(2))*edge(1,l) + &
                        si(2,face(2))*edge(2,l) + &
                        si(3,face(2))*edge(3,l))/(smag2*emag)
        ENDDO

        IF (MINVAL(dprod) < dprodm) THEN
          dprodm = MINVAL( dprod )
          im = i
          jm = j
          km = k
        ENDIF

        IF (MAXVAL(dprod) > dpmax) THEN
          dpmax = MAXVAL( dprod )
        ENDIF

! ----- check for skewed j-faces

        edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(3))
        edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(3))
        edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(3))

        edge(1,2) = xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(8))
        edge(2,2) = xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(8))
        edge(3,2) = xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(8))

        edge(1,3) = xyz(XCOORD,ijkNode(7))-xyz(XCOORD,ijkNode(5))
        edge(2,3) = xyz(YCOORD,ijkNode(7))-xyz(YCOORD,ijkNode(5))
        edge(3,3) = xyz(ZCOORD,ijkNode(7))-xyz(ZCOORD,ijkNode(5))

        edge(1,4) = xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(6))
        edge(2,4) = xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(6))
        edge(3,4) = xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(6))

        DO l=1,4
          emag  = SQRT( edge(1,l)*edge(1,l) + &
                        edge(2,l)*edge(2,l) + &
                        edge(3,l)*edge(3,l) )
          smag1 = SQRT( sj(1,face(1))*sj(1,face(1)) + &
                        sj(2,face(1))*sj(2,face(1)) + &
                        sj(3,face(1))*sj(3,face(1)) )
          smag2 = SQRT( sj(1,face(4))*sj(1,face(4)) + &
                        sj(2,face(4))*sj(2,face(4)) + &
                        sj(3,face(4))*sj(3,face(4)) )

          dprod(l)  =  (sj(1,face(1))*edge(1,l) + &
                        sj(2,face(1))*edge(2,l) + &
                        sj(3,face(1))*edge(3,l))/(smag1*emag)

          dprod(l+4) = (sj(1,face(4))*edge(1,l) + &
                        sj(2,face(4))*edge(2,l) + &
                        sj(3,face(4))*edge(3,l))/(smag2*emag)
        ENDDO

        IF (MINVAL(dprod) < dprodm) THEN
          dprodm = MINVAL( dprod )
          im = i
          jm = j
          km = k
        ENDIF

        IF (MAXVAL(dprod) > dpmax) THEN
          dpmax = MAXVAL( dprod )
        ENDIF

! ----- check for skewed k-faces

        edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(4))
        edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(4))
        edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(4))

        edge(1,2) = xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(7))
        edge(2,2) = xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(7))
        edge(3,2) = xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(7))

        edge(1,3) = xyz(XCOORD,ijkNode(6))-xyz(XCOORD,ijkNode(5))
        edge(2,3) = xyz(YCOORD,ijkNode(6))-xyz(YCOORD,ijkNode(5))
        edge(3,3) = xyz(ZCOORD,ijkNode(6))-xyz(ZCOORD,ijkNode(5))

        edge(1,4) = xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(8))
        edge(2,4) = xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(8))
        edge(3,4) = xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(8))

        DO l=1,4
          emag  = SQRT( edge(1,l)*edge(1,l) + &
                        edge(2,l)*edge(2,l) + &
                        edge(3,l)*edge(3,l) )
          smag1 = SQRT( sk(1,face(1))*sk(1,face(1)) + &
                        sk(2,face(1))*sk(2,face(1)) + &
                        sk(3,face(1))*sk(3,face(1)) )
          smag2 = SQRT( sk(1,face(6))*sk(1,face(6)) + &
                        sk(2,face(6))*sk(2,face(6)) + &
                        sk(3,face(6))*sk(3,face(6)) )

          dprod(l)  =  (sk(1,face(1))*edge(1,l) + &
                        sk(2,face(1))*edge(2,l) + &
                        sk(3,face(1))*edge(3,l))/(smag1*emag)

          dprod(l+4) = (sk(1,face(6))*edge(1,l) + &
                        sk(2,face(6))*edge(2,l) + &
                        sk(3,face(6))*edge(3,l))/(smag2*emag)
        ENDDO

        IF (MINVAL(dprod) < dprodm) THEN
          dprodm = MINVAL( dprod )
          im = i
          jm = j
          km = k
        ENDIF

        IF (MAXVAL(dprod) > dpmax) THEN
          dpmax = MAXVAL( dprod )
        ENDIF

      ENDDO
    ENDDO
  ENDDO

  IF (dprodm < 1.E-1_RFREAL) THEN
    region%blockShape = REGION_SHAPE_FUNKY

    IF ( region%global%verbLevel >= VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,I6)') SOLVER_NAME, &
           'Warning: found funky block, region',region%iRegionGlobal
    END IF ! global%verbLevel
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_BlockSkewedCell


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModGridRegionShape


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModGridRegionShape.F90,v $
! Revision 1.3  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/03/18 11:09:30  wasistho
! initial import
!
!
!
! ******************************************************************************









