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
! Purpose: Suite for grid smoothing routines based on elliptic PDE.
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModElliptSmoothing.F90,v 1.15 2008/12/06 08:44:15 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModElliptSmoothing

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
  PUBLIC :: RFLO_ElliptGridSmoo, &
            RFLO_ElliptGridSmooRegion, &
            RFLO_ElliptGridPatch, &
            RFLO_ElliptGridJac3D, &
            RFLO_ElliptGridJac2D, &
            RFLO_ElliptGridGauss3D, &
            RFLO_ElliptGridGauss2D, &
            RFLO_ElliptGridSOR3D, &
            RFLO_ElliptGridSOR2D, &
            RFLO_ElliptMatrixCoeffs3D, &
            RFLO_ElliptMatrixCoeffs2D

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModElliptSmoothing.F90,v $ $Revision: 1.15 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: smooth the distribution of grid points by solving 3D elliptic PDE 
!          in physical space.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%grid%xyz = new grid coordinates
!         resid = convergence of the Jacobi iteration.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridSmoo( regions,resid )

  USE ModInterfaces, ONLY : RFLO_GenerateCoarseGrids, &
        RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry, &
        RFLO_ExchangeGeometry, RFLO_ExchangeDnodeCopy, &
        RFLO_ExchangeDnodeSend, RFLO_ExchangeDnodeRecv, RFLO_ClearSendRequests,&
        RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset, RFLO_CalcFaceVectors, &
        RFLO_CalcCellCentroids, RFLO_CalcFaceCentroids, RFLO_ChangeInteriorGrid
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridPhysGrad3D
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  REAL(RFREAL) :: resid
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ijk, i, j, k

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: bcType, iRegSrc, iPatchSrc

  REAL(RFREAL) :: dx, dy, dz
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:), xyzOrig(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER   :: grid, gridOld, gridSrc
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ElliptGridSmoo',&
       'RFLO_ModElliptSmoothing.F90' )

! smooth grid region-wise -----------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      CALL RFLO_GridPhysGrad3D( regions(iReg) )

! --- copy coordinates

      xyz    => regions(iReg)%levels(1)%grid%xyz
      xyzOld => regions(iReg)%levels(1)%grid%xyzOld
      xyzOld =  xyz

! --- compute grid control map

      CALL RFLO_ElliptMatrixCoeffs3D( regions(iReg) )

      CALL RFLO_ElliptGridJac3D( regions(iReg) )
!      CALL RFLO_ElliptGridGauss3D( regions(iReg) )
!      CALL RFLO_ElliptGridSOR3D( regions(iReg) )

! --- transform grid coordinates to deformations

      xyz     => regions(iReg)%levels(1)%grid%xyz
      xyzOrig => regions(iReg)%levels(1)%gridOld%xyz

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
      ENDDO

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyzOld(XCOORD,ijk) = xyzOld(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
        xyzOld(YCOORD,ijk) = xyzOld(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
        xyzOld(ZCOORD,ijk) = xyzOld(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
      ENDDO

! --- zero out movements along certain boundaries, comment it out
!     in case deformations extended to boundary points

!      DO iPatch=1,regions(iReg)%nPatches
!        patch  => regions(iReg)%levels(1)%patches(iPatch)
!        bcType =  patch%bcType
!        IF ((bcType>=BC_INFLOW     .AND. bcType<=BC_INFLOW    +BC_RANGE) .OR. &
!            (bcType>=BC_OUTFLOW    .AND. bcType<=BC_OUTFLOW   +BC_RANGE) .OR. &
!            (bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL  +BC_RANGE) .OR. &
!            (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
!            (bcType>=BC_FARFIELD   .AND. bcType<=BC_FARFIELD  +BC_RANGE) .OR. &
!            (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION +BC_RANGE) .OR. &
!            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
!            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
!          CALL RFLO_ElliptGridPatch( regions(iReg),patch )
!        ENDIF  ! bcType
!      ENDDO    ! iPatch

      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      CALL RFLO_GetNodeOffset( regions(iReg),1,iNOff,ijNOff )

      resid = 0._RFREAL
      DO k=kpnbeg,kpnend
        DO j=jpnbeg,jpnend
          DO i=ipnbeg,ipnend
            ijk   = IndIJK(i,j,k,iNOff,ijNOff)
            dx    = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
            dy    = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
            dz    = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)

            resid = resid + dx*dx + dy*dy +dz*dz
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! fix interfaces between regions ----------------------------------------------
! copy / send deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid == global%myProcid) THEN
            CALL RFLO_ExchangeDnodeCopy( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true., &
                                         grid%xyz,gridSrc%xyz )
          ELSE
            CALL RFLO_ExchangeDnodeSend( regions(iReg),regions(iRegSrc), &
                                         patch,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! receive deformations

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

      grid    => regions(iReg)%levels(1)%grid
      gridOld => regions(iReg)%levels(1)%gridOld

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(1)%patches(iPatch)
        bcType =  patch%bcType
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
          iRegSrc   =  patch%srcRegion
          iPatchSrc =  patch%srcPatch
          patchSrc  => regions(iRegSrc)%levels(1)%patches(iPatchSrc)
          gridSrc   => regions(iRegSrc)%levels(1)%grid

          IF (regions(iRegSrc)%procid /= global%myProcid) THEN
            CALL RFLO_ExchangeDnodeRecv( regions(iReg),regions(iRegSrc), &
                                         patch,patchSrc,.true.,grid%xyz )
          ENDIF
        ENDIF  ! bcType
      ENDDO    ! iPatch

    ENDIF  ! region on this processor and active, grid moving
  ENDDO    ! iReg

! clear send requests

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving
      CALL RFLO_ClearSendRequests( regions,iReg,.true. )
    ENDIF
  ENDDO

! update grid, dummy, corner and edge cells -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE .AND. &           ! on my processor
        regions(iReg)%mixtInput%moveGrid) THEN         ! and moving

! --- change xyz from deformations to coordinates

      xyz     => regions(iReg)%levels(1)%grid%xyz
      xyzOrig => regions(iReg)%levels(1)%gridOld%xyz

      DO ijk=LBOUND(xyz,2),UBOUND(xyz,2)
        xyz(XCOORD,ijk) = xyz(XCOORD,ijk) + xyzOrig(XCOORD,ijk)
        xyz(YCOORD,ijk) = xyz(YCOORD,ijk) + xyzOrig(YCOORD,ijk)
        xyz(ZCOORD,ijk) = xyz(ZCOORD,ijk) + xyzOrig(ZCOORD,ijk)
      ENDDO

! --- update coarse grids and dummy cells

      CALL RFLO_GenerateCoarseGrids( regions(iReg) )   ! coarsen finest grid
      CALL RFLO_CopyGeometryDummy( regions(iReg) )     ! copy to dummy nodes
      CALL RFLO_ExtrapolateGeometry( regions(iReg) )   ! extrapolate
    ENDIF     ! region on this processor and active, grid moving
  ENDDO       ! iReg

  CALL RFLO_ExchangeGeometry( regions )                ! exchange geometry

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ElliptGridSmoo


!******************************************************************************
!
! Purpose: smooth the distribution of grid points regionwise by solving 
!          3D elliptic PDE in physical space.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%grid%xyz = new grid coordinates
!         resid = convergence of the Jacobi iteration.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridSmooRegion( region,resid )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridPhysGrad3D

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  REAL(RFREAL)   :: resid
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijk, i, j, k

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff

  REAL(RFREAL) :: dx, dy, dz
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridSmooRegion',&
       'RFLO_ModElliptSmoothing.F90' )

! compute gradients of grid coordinates ---------------------------------------

  CALL RFLO_GridPhysGrad3D( region )

! copy coordinates

  xyz    => region%levels(1)%grid%xyz
  xyzOld => region%levels(1)%grid%xyzOld
  xyzOld =  xyz

! compute grid control map

  CALL RFLO_ElliptMatrixCoeffs3D( region )

! solve for xyz

  CALL RFLO_ElliptGridJac3D( region )
!  CALL RFLO_ElliptGridGauss3D( region )
!  CALL RFLO_ElliptGridSOR3D( region )

! compute residual

  CALL RFLO_GetDimensPhysNodes( region,1,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,1,iNOff,ijNOff )

  resid = 0._RFREAL
  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk   = IndIJK(i,j,k,iNOff,ijNOff)
        dx    = xyz(XCOORD,ijk) - xyzOld(XCOORD,ijk)
        dy    = xyz(YCOORD,ijk) - xyzOld(YCOORD,ijk)
        dz    = xyz(ZCOORD,ijk) - xyzOld(ZCOORD,ijk)

        resid = resid + dx*dx + dy*dy +dz*dz
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ElliptGridSmooRegion



!******************************************************************************
!
! Purpose: zero out movements obtained by previous treatments along a boundary.
!
! Description: none.
!
! Input: region = data of current region, grid movements
!        patch  = current patch.
!
! Output: region%levels%grid%xyz = new grid movements.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridPatch( region,patch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff, ijkNB

  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ElliptGridPatch',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers

  iLev = 1

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                  jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz    => region%levels(iLev)%grid%xyz
  xyzOld => region%levels(iLev)%grid%xyzOld

! for grid coordinates: new = old

!  DO k=kbeg,kend
!    DO j=jbeg,jend
!      DO i=ibeg,iend
!        ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
!        xyz(XCOORD,ijkNB) = xyzOld(XCOORD,ijkNB)
!        xyz(YCOORD,ijkNB) = xyzOld(YCOORD,ijkNB)
!        xyz(ZCOORD,ijkNB) = xyzOld(ZCOORD,ijkNB)
!      ENDDO
!    ENDDO
!  ENDDO

! for deformations: zero out deformations

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkNB  = IndIJK(i,j,k,iNOff,ijNOff)
        xyz(XCOORD,ijkNB) = 0._RFREAL
        xyz(YCOORD,ijkNB) = 0._RFREAL
        xyz(ZCOORD,ijkNB) = 0._RFREAL
      ENDDO
    ENDDO
  ENDDO

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridPatch


!******************************************************************************
!
! Purpose: compute matrix coefficients for elliptic PDE grid smooting.
!
! Description: Equation 4.68 in Handbook of Grid Generation, by
!              Joe F. Thompson, Bharat K. Soni, and Nigel P. Weatherill,
!              CRC Press 1999, ISBN 0-8493-2687-7. 
!
! Input: region = grid data of current region
!
! Output: grid%aijk, aimjk, ... = matrix coefficients
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptMatrixCoeffs3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                            RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: iLev, ijkN, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend

  REAL(RFREAL) :: di, dj, dk, rdi2, rdj2, rdk2, r2di, r2dj, r2dk, q1, q2, q3
  REAL(RFREAL) :: al11, al12, al13, al22, al23, al33
  REAL(RFREAL) :: au11, au12, au13, au22, au23, au33
  REAL(RFREAL), POINTER :: stui(:,:), stuj(:,:), stuk(:,:), pmat(:,:,:)
  REAL(RFREAL), POINTER :: aijk(:), aimjk(:), aipjk(:), aijmk(:), aijpk(:)
  REAL(RFREAL), POINTER :: aijkm(:), aijkp(:), aimjmk(:), aipjmk(:)
  REAL(RFREAL), POINTER :: aimjpk(:), aipjpk(:), aimjkm(:), aipjkm(:)
  REAL(RFREAL), POINTER :: aimjkp(:), aipjkp(:), aijmkm(:), aijpkm(:)
  REAL(RFREAL), POINTER :: aijmkp(:), aijpkp(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptMatrixCoeffs3D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  stui   => region%levels(iLev)%grid%stui
  stuj   => region%levels(iLev)%grid%stuj
  stuk   => region%levels(iLev)%grid%stuk
  pmat   => region%levels(iLev)%grid%pmat

  aijk   => region%levels(iLev)%grid%aijk   
  aimjk  => region%levels(iLev)%grid%aimjk  
  aipjk  => region%levels(iLev)%grid%aipjk  
  aijmk  => region%levels(iLev)%grid%aijmk  
  aijpk  => region%levels(iLev)%grid%aijpk  
  aijkm  => region%levels(iLev)%grid%aijkm  
  aijkp  => region%levels(iLev)%grid%aijkp  
  aimjmk => region%levels(iLev)%grid%aimjmk 
  aipjmk => region%levels(iLev)%grid%aipjmk 
  aimjpk => region%levels(iLev)%grid%aimjpk 
  aipjpk => region%levels(iLev)%grid%aipjpk 
  aimjkm => region%levels(iLev)%grid%aimjkm 
  aipjkm => region%levels(iLev)%grid%aipjkm 
  aimjkp => region%levels(iLev)%grid%aimjkp 
  aipjkp => region%levels(iLev)%grid%aipjkp 
  aijmkm => region%levels(iLev)%grid%aijmkm 
  aijpkm => region%levels(iLev)%grid%aijpkm 
  aijmkp => region%levels(iLev)%grid%aijmkp 
  aijpkp => region%levels(iLev)%grid%aijpkp 

  di = 1._RFREAL/(ipnend-ipnbeg)
  dj = 1._RFREAL/(jpnend-jpnbeg)
  dk = 1._RFREAL/(kpnend-kpnbeg)

  rdi2 = 1._RFREAL/(di*di)
  rdj2 = 1._RFREAL/(dj*dj)
  rdk2 = 1._RFREAL/(dk*dk)
  r2di = 1._RFREAL/(2*di)
  r2dj = 1._RFREAL/(2*dj)
  r2dk = 1._RFREAL/(2*dk)

! start -----------------------------------------------------------------------

  DO k=kdnbeg,kdnend
    DO j=jdnbeg,jdnend
      DO i=idnbeg,idnend
        ijkN   = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)

        al11 = DOT_PRODUCT( stui(:,ijkN),stui(:,ijkN) ) 
        al12 = DOT_PRODUCT( stui(:,ijkN),stuj(:,ijkN) ) 
        al13 = DOT_PRODUCT( stui(:,ijkN),stuk(:,ijkN) ) 
        al22 = DOT_PRODUCT( stuj(:,ijkN),stuj(:,ijkN) ) 
        al23 = DOT_PRODUCT( stuj(:,ijkN),stuk(:,ijkN) ) 
        al33 = DOT_PRODUCT( stuk(:,ijkN),stuk(:,ijkN) ) 
        au11 = al22*al33 - al23*al23
        au12 = al13*al23 - al12*al33
        au13 = al12*al23 - al13*al22
        au22 = al11*al33 - al13*al13
        au23 = al13*al12 - al11*al23
        au33 = al11*al22 - al12*al12

        q1 = au11*pmat(XCOORD,1,ijkN) + 2._RFREAL*au12*pmat(XCOORD,2,ijkN) + &
             2._RFREAL*au13*pmat(XCOORD,3,ijkN) + au22*pmat(XCOORD,4,ijkN) + &
             2._RFREAL*au23*pmat(XCOORD,5,ijkN) + au33*pmat(XCOORD,6,ijkN)
        q2 = au11*pmat(YCOORD,1,ijkN) + 2._RFREAL*au12*pmat(YCOORD,2,ijkN) + &
             2._RFREAL*au13*pmat(YCOORD,3,ijkN) + au22*pmat(YCOORD,4,ijkN) + &
             2._RFREAL*au23*pmat(YCOORD,5,ijkN) + au33*pmat(YCOORD,6,ijkN)
        q3 = au11*pmat(ZCOORD,1,ijkN) + 2._RFREAL*au12*pmat(ZCOORD,2,ijkN) + &
             2._RFREAL*au13*pmat(ZCOORD,3,ijkN) + au22*pmat(ZCOORD,4,ijkN) + &
             2._RFREAL*au23*pmat(ZCOORD,5,ijkN) + au33*pmat(ZCOORD,6,ijkN)

        aijk(  ijkN) =  2._RFREAL*(au11*rdi2 + au22*rdj2 + au33*rdk2)
        aimjk( ijkN) =  au11*rdi2 - q1*r2di
        aipjk( ijkN) =  au11*rdi2 + q1*r2di
        aijmk( ijkN) =  au22*rdj2 - q2*r2dj 
        aijpk( ijkN) =  au22*rdj2 + q2*r2dj 
        aijkm( ijkN) =  au33*rdk2 - q3*r2dk  
        aijkp( ijkN) =  au33*rdk2 + q3*r2dk 
        aimjmk(ijkN) =  2._RFREAL*au12*r2di*r2dj 
        aipjmk(ijkN) = -2._RFREAL*au12*r2di*r2dj 
        aimjpk(ijkN) = -2._RFREAL*au12*r2di*r2dj  
        aipjpk(ijkN) =  2._RFREAL*au12*r2di*r2dj  
        aimjkm(ijkN) =  2._RFREAL*au13*r2di*r2dk  
        aipjkm(ijkN) = -2._RFREAL*au13*r2di*r2dk  
        aimjkp(ijkN) = -2._RFREAL*au13*r2di*r2dk  
        aipjkp(ijkN) =  2._RFREAL*au13*r2di*r2dk  
        aijmkm(ijkN) =  2._RFREAL*au23*r2dj*r2dk  
        aijpkm(ijkN) = -2._RFREAL*au23*r2dj*r2dk  
        aijmkp(ijkN) = -2._RFREAL*au23*r2dj*r2dk  
        aijpkp(ijkN) =  2._RFREAL*au23*r2dj*r2dk  
!$BTEST
!        write(*,*) 'pmat1',region%iRegionGlobal,i,j,k,pmat(1:3,1,ijkN)
!        write(*,*) 'pmat2',region%iRegionGlobal,i,j,k,pmat(1:3,2,ijkN)
!        write(*,*) 'pmat3',region%iRegionGlobal,i,j,k,pmat(1:3,3,ijkN)
!        write(*,*) 'pmat4',region%iRegionGlobal,i,j,k,pmat(1:3,4,ijkN)
!        write(*,*) 'pmat5',region%iRegionGlobal,i,j,k,pmat(1:3,5,ijkN)
!        write(*,*) 'pmat6',region%iRegionGlobal,i,j,k,pmat(1:3,6,ijkN)
!$ETEST
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptMatrixCoeffs3D


!******************************************************************************
!
! Purpose: compute matrix coefficients for elliptic PDE grid smooting.
!
! Description: Equation 4.20 in Handbook of Grid Generation, by
!              Joe F. Thompson, Bharat K. Soni, and Nigel P. Weatherill,
!              CRC Press 1999, ISBN 0-8493-2687-7. 
!
! Input: region = grid data of current region
!
! Output: grid%aijk, aimjk, ... = matrix coefficients
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptMatrixCoeffs2D( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
                            RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: h1, h2

  REAL(RFREAL) :: di, dj, rdi2, rdj2, r2di, r2dj, pco, qco, rco, sco, tco
  REAL(RFREAL), POINTER :: sti(:,:,:), stj(:,:,:), pfun(:,:,:,:)
  REAL(RFREAL), POINTER :: aij(:,:), aimj(:,:), aipj(:,:), aijm(:,:), aijp(:,:)
  REAL(RFREAL), POINTER :: aimjm(:,:), aipjm(:,:), aimjp(:,:), aipjp(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptMatrixCoeffs2D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  sti   => patch%sti
  stj   => patch%stj
  pfun  => patch%pfun

  aij   => patch%aij  
  aimj  => patch%aimj 
  aipj  => patch%aipj 
  aijm  => patch%aijm 
  aijp  => patch%aijp 
  aimjm => patch%aimjm
  aipjm => patch%aipjm
  aimjp => patch%aimjp
  aipjp => patch%aipjp

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2
  
  di = 1._RFREAL/(h1-1)
  dj = 1._RFREAL/(h2-1)

  rdi2 = 1._RFREAL/(di*di)
  rdj2 = 1._RFREAL/(dj*dj)
  r2di = 1._RFREAL/(2*di)
  r2dj = 1._RFREAL/(2*dj)

! start -----------------------------------------------------------------------

  DO j=1,h2
    DO i=1,h1
      pco = DOT_PRODUCT( stj(:,i,j),stj(:,i,j) ) 
      qco = DOT_PRODUCT( sti(:,i,j),stj(:,i,j) ) 
      rco = DOT_PRODUCT( sti(:,i,j),sti(:,i,j) ) 
      sco = pco*pfun(1,1,i,j) - 2._RFREAL*qco*pfun(1,2,i,j) + rco*pfun(1,3,i,j)
      tco = pco*pfun(2,1,i,j) - 2._RFREAL*qco*pfun(2,2,i,j) + rco*pfun(2,3,i,j)

      aij(  i,j) =   2._RFREAL*(pco*rdi2 + rco*rdj2)
      aimj( i,j) =   pco*rdi2 - sco*r2di
      aipj( i,j) =   pco*rdi2 + sco*r2di
      aijm( i,j) =   rco*rdj2 - tco*r2dj
      aijp( i,j) =   rco*rdj2 + tco*r2dj
      aimjm(i,j) =  -2._RFREAL*qco*r2di*r2dj 
      aipjm(i,j) =   2._RFREAL*qco*r2di*r2dj 
      aimjp(i,j) =   2._RFREAL*qco*r2di*r2dj 
      aipjp(i,j) =  -2._RFREAL*qco*r2di*r2dj 
!$BTEST
!      IF ((region%iRegionGlobal == 24 .OR. &
!           region%iRegionGlobal == 27 .OR. &
!           region%iRegionGlobal == 30 .OR. &
!           region%iRegionGlobal == 33) .AND. iPatch==2) &
!        write(*,*)region%iRegionGlobal,iPatch,i,j, &
!        pfun(1:2,1,i,j), &
!        pfun(1:2,2,i,j), &
!        pfun(1:2,3,i,j)
!$ETEST
    ENDDO ! i
  ENDDO   ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptMatrixCoeffs2D


!******************************************************************************
!
! Purpose: perform one Jacobi iteration to solve for physical space grid
!          coordinate by elliptic PDE.
!
! Description: the method used is Gauss-Jacobi method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridJac3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijkN, imjkN, ipjkN, ijmkN, ijpkN, ijkmN, ijkpN
  INTEGER :: imjmkN, ipjmkN, imjpkN, ipjpkN, imjkmN, ipjkmN, imjkpN, ipjkpN
  INTEGER :: ijmkmN, ijpkmN, ijmkpN, ijpkpN

  REAL(RFREAL) :: raijk
  REAL(RFREAL), POINTER :: aijk(:), aimjk(:), aipjk(:), aijmk(:), aijpk(:)
  REAL(RFREAL), POINTER :: aijkm(:), aijkp(:), aimjmk(:), aipjmk(:), aimjpk(:)
  REAL(RFREAL), POINTER :: aipjpk(:), aimjkm(:), aipjkm(:), aimjkp(:)
  REAL(RFREAL), POINTER :: aipjkp(:), aijmkm(:), aijpkm(:), aijmkp(:)
  REAL(RFREAL), POINTER :: aijpkp(:), xyz(:,:), xyzOld(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridJac3D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  aijk   => region%levels(iLev)%grid%aijk   
  aimjk  => region%levels(iLev)%grid%aimjk  
  aipjk  => region%levels(iLev)%grid%aipjk  
  aijmk  => region%levels(iLev)%grid%aijmk  
  aijpk  => region%levels(iLev)%grid%aijpk  
  aijkm  => region%levels(iLev)%grid%aijkm  
  aijkp  => region%levels(iLev)%grid%aijkp  
  aimjmk => region%levels(iLev)%grid%aimjmk 
  aipjmk => region%levels(iLev)%grid%aipjmk 
  aimjpk => region%levels(iLev)%grid%aimjpk 
  aipjpk => region%levels(iLev)%grid%aipjpk 
  aimjkm => region%levels(iLev)%grid%aimjkm 
  aipjkm => region%levels(iLev)%grid%aipjkm 
  aimjkp => region%levels(iLev)%grid%aimjkp 
  aipjkp => region%levels(iLev)%grid%aipjkp 
  aijmkm => region%levels(iLev)%grid%aijmkm 
  aijpkm => region%levels(iLev)%grid%aijpkm 
  aijmkp => region%levels(iLev)%grid%aijmkp 
  aijpkp => region%levels(iLev)%grid%aijpkp 

  xyz    => region%levels(iLev)%grid%xyz
  xyzOld => region%levels(iLev)%grid%xyzOld

! perform Gauss-Jacobi --------------------------------------------------------

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ipjkN   = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijpkN   = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkpN   = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        imjmkN  = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)
        ipjmkN  = IndIJK(i+1,j-1,k  ,iNOff,ijNOff)
        imjpkN  = IndIJK(i-1,j+1,k  ,iNOff,ijNOff)
        ipjpkN  = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
        imjkmN  = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)
        ipjkmN  = IndIJK(i+1,j  ,k-1,iNOff,ijNOff)
        imjkpN  = IndIJK(i-1,j  ,k+1,iNOff,ijNOff)
        ipjkpN  = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        ijmkmN  = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)
        ijpkmN  = IndIJK(i  ,j+1,k-1,iNOff,ijNOff)
        ijmkpN  = IndIJK(i  ,j-1,k+1,iNOff,ijNOff)
        ijpkpN  = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

        raijk   = 1._RFREAL/aijk(ijkN)

        DO l = XCOORD,ZCOORD
          xyz(l,ijkN) = aimjk(  ijkN)*xyzOld(l, imjkN) + &
                        aipjk(  ijkN)*xyzOld(l, ipjkN) + &
                        aijmk(  ijkN)*xyzOld(l, ijmkN) + &
                        aijpk(  ijkN)*xyzOld(l, ijpkN) + &
                        aijkm(  ijkN)*xyzOld(l, ijkmN) + &
                        aijkp(  ijkN)*xyzOld(l, ijkpN) + &
                        aimjmk( ijkN)*xyzOld(l,imjmkN) + &
                        aipjmk( ijkN)*xyzOld(l,ipjmkN) + &
                        aimjpk( ijkN)*xyzOld(l,imjpkN) + &
                        aipjpk( ijkN)*xyzOld(l,ipjpkN) + &
                        aimjkm( ijkN)*xyzOld(l,imjkmN) + &
                        aipjkm( ijkN)*xyzOld(l,ipjkmN) + &
                        aimjkp( ijkN)*xyzOld(l,imjkpN) + &
                        aipjkp( ijkN)*xyzOld(l,ipjkpN) + &
                        aijmkm( ijkN)*xyzOld(l,ijmkmN) + &
                        aijpkm( ijkN)*xyzOld(l,ijpkmN) + &
                        aijmkp( ijkN)*xyzOld(l,ijmkpN) + &
                        aijpkp( ijkN)*xyzOld(l,ijpkpN)
          xyz(l,ijkN) = xyz(l,ijkN)*raijk
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridJac3D


!******************************************************************************
!
! Purpose: perform one Gauss-Seidel iteration to solve for physical space grid
!          coordinate by elliptic PDE.
!
! Description: the method used is Gauss-Seidel method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridGauss3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijkN, imjkN, ipjkN, ijmkN, ijpkN, ijkmN, ijkpN
  INTEGER :: imjmkN, ipjmkN, imjpkN, ipjpkN, imjkmN, ipjkmN, imjkpN, ipjkpN
  INTEGER :: ijmkmN, ijpkmN, ijmkpN, ijpkpN

  REAL(RFREAL) :: raijk
  REAL(RFREAL), POINTER :: aijk(:), aimjk(:), aipjk(:), aijmk(:), aijpk(:)
  REAL(RFREAL), POINTER :: aijkm(:), aijkp(:), aimjmk(:), aipjmk(:), aimjpk(:)
  REAL(RFREAL), POINTER :: aipjpk(:), aimjkm(:), aipjkm(:), aimjkp(:)
  REAL(RFREAL), POINTER :: aipjkp(:), aijmkm(:), aijpkm(:), aijmkp(:)
  REAL(RFREAL), POINTER :: aijpkp(:), xyz(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridGauss3D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  aijk   => region%levels(iLev)%grid%aijk   
  aimjk  => region%levels(iLev)%grid%aimjk  
  aipjk  => region%levels(iLev)%grid%aipjk  
  aijmk  => region%levels(iLev)%grid%aijmk  
  aijpk  => region%levels(iLev)%grid%aijpk  
  aijkm  => region%levels(iLev)%grid%aijkm  
  aijkp  => region%levels(iLev)%grid%aijkp  
  aimjmk => region%levels(iLev)%grid%aimjmk 
  aipjmk => region%levels(iLev)%grid%aipjmk 
  aimjpk => region%levels(iLev)%grid%aimjpk 
  aipjpk => region%levels(iLev)%grid%aipjpk 
  aimjkm => region%levels(iLev)%grid%aimjkm 
  aipjkm => region%levels(iLev)%grid%aipjkm 
  aimjkp => region%levels(iLev)%grid%aimjkp 
  aipjkp => region%levels(iLev)%grid%aipjkp 
  aijmkm => region%levels(iLev)%grid%aijmkm 
  aijpkm => region%levels(iLev)%grid%aijpkm 
  aijmkp => region%levels(iLev)%grid%aijmkp 
  aijpkp => region%levels(iLev)%grid%aijpkp 

  xyz    => region%levels(iLev)%grid%xyz

! perform Gauss iteration ----------------------------------------------------

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ipjkN   = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijpkN   = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkpN   = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        imjmkN  = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)
        ipjmkN  = IndIJK(i+1,j-1,k  ,iNOff,ijNOff)
        imjpkN  = IndIJK(i-1,j+1,k  ,iNOff,ijNOff)
        ipjpkN  = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
        imjkmN  = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)
        ipjkmN  = IndIJK(i+1,j  ,k-1,iNOff,ijNOff)
        imjkpN  = IndIJK(i-1,j  ,k+1,iNOff,ijNOff)
        ipjkpN  = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        ijmkmN  = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)
        ijpkmN  = IndIJK(i  ,j+1,k-1,iNOff,ijNOff)
        ijmkpN  = IndIJK(i  ,j-1,k+1,iNOff,ijNOff)
        ijpkpN  = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

        raijk   = 1._RFREAL/aijk(ijkN)

        DO l = XCOORD,ZCOORD
          xyz(l,ijkN) = aimjk(  ijkN)*xyz(l, imjkN) + &
                        aipjk(  ijkN)*xyz(l, ipjkN) + &
                        aijmk(  ijkN)*xyz(l, ijmkN) + &
                        aijpk(  ijkN)*xyz(l, ijpkN) + &
                        aijkm(  ijkN)*xyz(l, ijkmN) + &
                        aijkp(  ijkN)*xyz(l, ijkpN) + &
                        aimjmk( ijkN)*xyz(l,imjmkN) + &
                        aipjmk( ijkN)*xyz(l,ipjmkN) + &
                        aimjpk( ijkN)*xyz(l,imjpkN) + &
                        aipjpk( ijkN)*xyz(l,ipjpkN) + &
                        aimjkm( ijkN)*xyz(l,imjkmN) + &
                        aipjkm( ijkN)*xyz(l,ipjkmN) + &
                        aimjkp( ijkN)*xyz(l,imjkpN) + &
                        aipjkp( ijkN)*xyz(l,ipjkpN) + &
                        aijmkm( ijkN)*xyz(l,ijmkmN) + &
                        aijpkm( ijkN)*xyz(l,ijpkmN) + &
                        aijmkp( ijkN)*xyz(l,ijmkpN) + &
                        aijpkp( ijkN)*xyz(l,ijpkpN)
          xyz(l,ijkN) = xyz(l,ijkN)*raijk
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridGauss3D


!******************************************************************************
!
! Purpose: perform one SOR iteration to solve for physical space grid
!          coordinate by elliptic PDE.
!
! Description: the method used is SOR method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridSOR3D( region )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ijkN, imjkN, ipjkN, ijmkN, ijpkN, ijkmN, ijkpN
  INTEGER :: imjmkN, ipjmkN, imjpkN, ipjpkN, imjkmN, ipjkmN, imjkpN, ipjkpN
  INTEGER :: ijmkmN, ijpkmN, ijmkpN, ijpkpN

  REAL(RFREAL) :: raijk, xyzTemp, omega, omomg
  REAL(RFREAL), POINTER :: aijk(:), aimjk(:), aipjk(:), aijmk(:), aijpk(:)
  REAL(RFREAL), POINTER :: aijkm(:), aijkp(:), aimjmk(:), aipjmk(:), aimjpk(:)
  REAL(RFREAL), POINTER :: aipjpk(:), aimjkm(:), aipjkm(:), aimjkp(:)
  REAL(RFREAL), POINTER :: aipjkp(:), aijmkm(:), aijpkm(:), aijmkp(:)
  REAL(RFREAL), POINTER :: aijpkp(:), xyz(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridSOR3D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  aijk   => region%levels(iLev)%grid%aijk   
  aimjk  => region%levels(iLev)%grid%aimjk  
  aipjk  => region%levels(iLev)%grid%aipjk  
  aijmk  => region%levels(iLev)%grid%aijmk  
  aijpk  => region%levels(iLev)%grid%aijpk  
  aijkm  => region%levels(iLev)%grid%aijkm  
  aijkp  => region%levels(iLev)%grid%aijkp  
  aimjmk => region%levels(iLev)%grid%aimjmk 
  aipjmk => region%levels(iLev)%grid%aipjmk 
  aimjpk => region%levels(iLev)%grid%aimjpk 
  aipjpk => region%levels(iLev)%grid%aipjpk 
  aimjkm => region%levels(iLev)%grid%aimjkm 
  aipjkm => region%levels(iLev)%grid%aipjkm 
  aimjkp => region%levels(iLev)%grid%aimjkp 
  aipjkp => region%levels(iLev)%grid%aipjkp 
  aijmkm => region%levels(iLev)%grid%aijmkm 
  aijpkm => region%levels(iLev)%grid%aijpkm 
  aijmkp => region%levels(iLev)%grid%aijmkp 
  aijpkp => region%levels(iLev)%grid%aijpkp 

  xyz    => region%levels(iLev)%grid%xyz

! perform point SOR -----------------------------------------------------------

  omega = 1.3_RFREAL
  omomg = 1._RFREAL - omega

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1
        ijkN    = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        imjkN   = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        ipjkN   = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ijmkN   = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        ijpkN   = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        ijkmN   = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        ijkpN   = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        imjmkN  = IndIJK(i-1,j-1,k  ,iNOff,ijNOff)
        ipjmkN  = IndIJK(i+1,j-1,k  ,iNOff,ijNOff)
        imjpkN  = IndIJK(i-1,j+1,k  ,iNOff,ijNOff)
        ipjpkN  = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
        imjkmN  = IndIJK(i-1,j  ,k-1,iNOff,ijNOff)
        ipjkmN  = IndIJK(i+1,j  ,k-1,iNOff,ijNOff)
        imjkpN  = IndIJK(i-1,j  ,k+1,iNOff,ijNOff)
        ipjkpN  = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        ijmkmN  = IndIJK(i  ,j-1,k-1,iNOff,ijNOff)
        ijpkmN  = IndIJK(i  ,j+1,k-1,iNOff,ijNOff)
        ijmkpN  = IndIJK(i  ,j-1,k+1,iNOff,ijNOff)
        ijpkpN  = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

        raijk   = 1._RFREAL/aijk(ijkN)

        DO l = XCOORD,ZCOORD
          xyzTemp     = aimjk(  ijkN)*xyz(l, imjkN) + &
                        aipjk(  ijkN)*xyz(l, ipjkN) + &
                        aijmk(  ijkN)*xyz(l, ijmkN) + &
                        aijpk(  ijkN)*xyz(l, ijpkN) + &
                        aijkm(  ijkN)*xyz(l, ijkmN) + &
                        aijkp(  ijkN)*xyz(l, ijkpN) + &
                        aimjmk( ijkN)*xyz(l,imjmkN) + &
                        aipjmk( ijkN)*xyz(l,ipjmkN) + &
                        aimjpk( ijkN)*xyz(l,imjpkN) + &
                        aipjpk( ijkN)*xyz(l,ipjpkN) + &
                        aimjkm( ijkN)*xyz(l,imjkmN) + &
                        aipjkm( ijkN)*xyz(l,ipjkmN) + &
                        aimjkp( ijkN)*xyz(l,imjkpN) + &
                        aipjkp( ijkN)*xyz(l,ipjkpN) + &
                        aijmkm( ijkN)*xyz(l,ijmkmN) + &
                        aijpkm( ijkN)*xyz(l,ijpkmN) + &
                        aijmkp( ijkN)*xyz(l,ijmkpN) + &
                        aijpkp( ijkN)*xyz(l,ijpkpN)
          xyzTemp     = xyzTemp*raijk
          xyz(l,ijkN) = omega*xyzTemp+omomg*xyz(l,ijkN)
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridSOR3D


!******************************************************************************
!
! Purpose: perform Jacobi iteration to solve for physical space grid
!          coordinate on patches by elliptic PDE.
!
! Description: the method used is Gauss-Jacobi method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridJac2D( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridPhysGrad2D

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: i, j, k, l, iter

! ... local variables
  INTEGER :: iLev, lbound, h1, h2, ng1, ng2, ijkN
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff

  REAL(RFREAL) :: raij, dx, dy, dz, resid
  REAL(RFREAL), POINTER :: aij(:,:), aimj(:,:), aipj(:,:), aijm(:,:), aijp(:,:)
  REAL(RFREAL), POINTER :: aimjm(:,:), aipjm(:,:), aimjp(:,:), aipjp(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:), st(:,:,:), stOld(:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridJac2D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions --------------------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  xyz    => region%levels(iLev)%grid%xyz
  st     => patch%st
  stOld  => patch%stOld

  aij   => patch%aij  
  aimj  => patch%aimj 
  aipj  => patch%aipj 
  aijm  => patch%aijm 
  aijp  => patch%aijp 
  aimjm => patch%aimjm
  aipjm => patch%aipjm
  aimjp => patch%aimjp
  aipjp => patch%aipjp

! copy xyz to stOld -----------------------------------------------------------

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
        stOld(1,ng1,ng2) = xyz(XCOORD,ijkN)
        stOld(2,ng1,ng2) = xyz(YCOORD,ijkN)
        stOld(3,ng1,ng2) = xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

! perform Gauss-Jacobi --------------------------------------------------------

  DO iter = 1,global%moveGridSiter

    CALL RFLO_GridPhysGrad2D( region,patch )
    CALL RFLO_ElliptMatrixCoeffs2D( region,patch,iPatch )

    resid = 0._RFREAL
    DO j=2,h2-1
      DO i=2,h1-1

        raij = 1._RFREAL/aij(i,j)

        DO l = 1,3
          st(l,i,j) = aimj( i,j)*stOld(l,i-1,j  ) + &
                      aipj( i,j)*stOld(l,i+1,j  ) + &
                      aijm( i,j)*stOld(l,i  ,j-1) + &
                      aijp( i,j)*stOld(l,i  ,j+1) + &
                      aimjm(i,j)*stOld(l,i-1,j-1) + &
                      aipjm(i,j)*stOld(l,i+1,j-1) + &
                      aimjp(i,j)*stOld(l,i-1,j+1) + &
                      aipjp(i,j)*stOld(l,i+1,j+1)
          st(l,i,j) = st(l,i,j)*raij
        ENDDO ! l

        dx    = st(1,i,j) - stOld(1,i,j)
        dy    = st(2,i,j) - stOld(2,i,j)
        dz    = st(3,i,j) - stOld(3,i,j)
        resid = resid + dx*dx + dy*dy + dz*dz
      ENDDO   ! i
    ENDDO     ! j

!$BTEST
!    write(*,*)region%iRegionGlobal,iPatch,iter,SQRT(resid)
!$ETEST

    stOld = st

! - copy st back to xyz -------------------------------------------------------

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
          xyz(XCOORD,ijkN)    = st(1,ng1,ng2)   
          xyz(YCOORD,ijkN)    = st(2,ng1,ng2)   
          xyz(ZCOORD,ijkN)    = st(3,ng1,ng2)   
        ENDDO ! i
      ENDDO   ! j
    ENDDO     ! k

  ENDDO       ! iter

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridJac2D


!******************************************************************************
!
! Purpose: perform Gauss-Seidel iteration to solve for physical space grid
!          coordinate on patches by elliptic PDE.
!
! Description: the method used is Gauss-Seidel method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridGauss2D( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridPhysGrad2D

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: i, j, k, l, iter

! ... local variables
  INTEGER :: iLev, lbound, h1, h2, ng1, ng2, ijkN
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff

  REAL(RFREAL) :: raij, dx, dy, dz, resid
  REAL(RFREAL), POINTER :: aij(:,:), aimj(:,:), aipj(:,:), aijm(:,:), aijp(:,:)
  REAL(RFREAL), POINTER :: aimjm(:,:), aipjm(:,:), aimjp(:,:), aipjp(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:), st(:,:,:), stOld(:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridGauss2D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions --------------------------------------------------------------

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  xyz   => region%levels(iLev)%grid%xyz
  st    => patch%st
  stOld => patch%stOld

  aij   => patch%aij  
  aimj  => patch%aimj 
  aipj  => patch%aipj 
  aijm  => patch%aijm 
  aijp  => patch%aijp 
  aimjm => patch%aimjm
  aipjm => patch%aipjm
  aimjp => patch%aimjp
  aipjp => patch%aipjp

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
        st(1,ng1,ng2)    = xyz(XCOORD,ijkN)
        st(2,ng1,ng2)    = xyz(YCOORD,ijkN)
        st(3,ng1,ng2)    = xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

! perform Gauss-Seidel --------------------------------------------------------

  DO iter = 1,global%moveGridSiter

    CALL RFLO_GridPhysGrad2D( region,patch )
    CALL RFLO_ElliptMatrixCoeffs2D( region,patch,iPatch )

    stOld = st

    resid = 0._RFREAL
    DO j=2,h2-1
      DO i=2,h1-1

        raij = 1._RFREAL/aij(i,j)

        DO l = 1,3
          st(l,i,j) = aimj( i,j)*st(l,i-1,j  ) + &
                      aipj( i,j)*st(l,i+1,j  ) + &
                      aijm( i,j)*st(l,i  ,j-1) + &
                      aijp( i,j)*st(l,i  ,j+1) + &
                      aimjm(i,j)*st(l,i-1,j-1) + &
                      aipjm(i,j)*st(l,i+1,j-1) + &
                      aimjp(i,j)*st(l,i-1,j+1) + &
                      aipjp(i,j)*st(l,i+1,j+1)
          st(l,i,j) = st(l,i,j)*raij
        ENDDO ! l

        dx    = st(1,i,j) - stOld(1,i,j)
        dy    = st(2,i,j) - stOld(2,i,j)
        dz    = st(3,i,j) - stOld(3,i,j)
        resid = resid + dx*dx + dy*dy + dz*dz
      ENDDO   ! i
    ENDDO     ! j

!$BTEST
!    write(*,*)region%iRegionGlobal,iPatch,iter,SQRT(resid)
!$ETEST

! - copy st back to xyz -------------------------------------------------------

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
          xyz(XCOORD,ijkN)    = st(1,ng1,ng2)   
          xyz(YCOORD,ijkN)    = st(2,ng1,ng2)   
          xyz(ZCOORD,ijkN)    = st(3,ng1,ng2)   
        ENDDO ! i
      ENDDO   ! j
    ENDDO     ! k

  ENDDO       ! iter

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridGauss2D


!******************************************************************************
!
! Purpose: perform SOR iteration to solve for physical space grid
!          coordinate on patches by elliptic PDE.
!
! Description: the method used is SOR method.
!
! Input: region = grid data of current region
!
! Output: xyz = new grid in physical space
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ElliptGridSOR2D( region,patch,iPatch )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE RFLO_ModGridControlMap, ONLY : RFLO_GridPhysGrad2D

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER :: iPatch

! ... loop variables
  INTEGER :: i, j, k, l, iter

! ... local variables
  INTEGER :: iLev, lbound, h1, h2, ng1, ng2, ijkN
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff

  REAL(RFREAL) :: raij, stemp, omega, omomg, dx, dy, dz, resid
  REAL(RFREAL), POINTER :: aij(:,:), aimj(:,:), aipj(:,:), aijm(:,:), aijp(:,:)
  REAL(RFREAL), POINTER :: aimjm(:,:), aipjm(:,:), aimjp(:,:), aipjp(:,:)
  REAL(RFREAL), POINTER :: xyz(:,:), st(:,:,:), stOld(:,:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_ElliptGridSOR2D',&
       'RFLO_ModElliptSmoothing.F90' )

! get dimensions and parameters -----------------------------------------------

  omega = 1.5_RFREAL
  omomg = 1._RFREAL - omega

  iLev   = 1
  lbound = patch%lbound

  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  h1 = patch%l1end-patch%l1beg+2
  h2 = patch%l2end-patch%l2beg+2

! get pointers ----------------------------------------------------------------

  xyz   => region%levels(iLev)%grid%xyz
  st    => patch%st
  stOld => patch%stOld

  aij   => patch%aij  
  aimj  => patch%aimj 
  aipj  => patch%aipj 
  aijm  => patch%aijm 
  aijp  => patch%aijp 
  aimjm => patch%aimjm
  aipjm => patch%aipjm
  aimjp => patch%aimjp
  aipjp => patch%aipjp

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
        st(1,ng1,ng2)    = xyz(XCOORD,ijkN)
        st(2,ng1,ng2)    = xyz(YCOORD,ijkN)
        st(3,ng1,ng2)    = xyz(ZCOORD,ijkN)
      ENDDO
    ENDDO
  ENDDO

! perform point SOR ----------------------------------------------------------

  DO iter = 1,global%moveGridSiter

    CALL RFLO_GridPhysGrad2D( region,patch )
    CALL RFLO_ElliptMatrixCoeffs2D( region,patch,iPatch )

    stOld = st

    resid = 0._RFREAL
    DO j=2,h2-1
      DO i=2,h1-1

        raij = 1._RFREAL/aij(i,j)

        DO l = 1,3
          stemp = aimj( i,j)*st(l,i-1,j  ) + &
                  aipj( i,j)*st(l,i+1,j  ) + &
                  aijm( i,j)*st(l,i  ,j-1) + &
                  aijp( i,j)*st(l,i  ,j+1) + &
                  aimjm(i,j)*st(l,i-1,j-1) + &
                  aipjm(i,j)*st(l,i+1,j-1) + &
                  aimjp(i,j)*st(l,i-1,j+1) + &
                  aipjp(i,j)*st(l,i+1,j+1)
          stemp = stemp*raij
          st(l,i,j) = omega*stemp+omomg*st(l,i,j)
        ENDDO ! l

        dx    = st(1,i,j) - stOld(1,i,j)
        dy    = st(2,i,j) - stOld(2,i,j)
        dz    = st(3,i,j) - stOld(3,i,j)
        resid = resid + dx*dx + dy*dy + dz*dz
      ENDDO   ! i
    ENDDO     ! j

!$BTEST
!    write(*,*)region%iRegionGlobal,iPatch,iter,SQRT(resid)
!$ETEST

! - copy st back to xyz -------------------------------------------------------

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
          xyz(XCOORD,ijkN)    = st(1,ng1,ng2)   
          xyz(YCOORD,ijkN)    = st(2,ng1,ng2)   
          xyz(ZCOORD,ijkN)    = st(3,ng1,ng2)   
        ENDDO ! i
      ENDDO   ! j
    ENDDO     ! k

  ENDDO       ! iter

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ElliptGridSOR2D


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModElliptSmoothing


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModElliptSmoothing.F90,v $
! Revision 1.15  2008/12/06 08:44:15  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/03/14 04:37:13  wasistho
! removed obsolete RFLO_ElliptFlatPatch
!
! Revision 1.12  2006/03/12 10:28:33  wasistho
! defined boundFlat in ElliptFlatPatch
!
! Revision 1.11  2006/03/08 06:39:35  wasistho
! added moveGridSiter
!
! Revision 1.10  2006/03/05 19:06:22  wasistho
! set computational space coordinates from initial grid
!
! Revision 1.9  2006/03/03 04:17:00  wasistho
! elliptic smoothings on coords instead of deformations
!
! Revision 1.8  2006/03/03 00:47:41  wasistho
! added elliptGridSmooRegion
!
! Revision 1.7  2006/03/02 00:23:54  wasistho
! prepared elliptic pde grid motion
!
! Revision 1.6  2006/02/11 03:42:14  wasistho
! added ModMoveGridElliptGlo/Fra
!
! Revision 1.5  2006/02/09 00:25:11  wasistho
! removed unused pointer xyzTemp
!
! Revision 1.4  2006/02/08 07:51:56  wasistho
! increased eps in flatPatch
!
! Revision 1.3  2005/12/07 08:46:15  wasistho
! added stuff for surface mesh motion EPDE
!
! Revision 1.2  2005/12/05 10:49:44  wasistho
! added RFLO_ElliptFlatPatch
!
! Revision 1.1  2005/12/03 09:39:47  wasistho
! initial import
!
!
!
! ******************************************************************************

















