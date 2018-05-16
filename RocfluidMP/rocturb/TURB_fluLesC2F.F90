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
! Purpose: Get face values of cv for unstructured grid.
!
! Description: It is obtained by two point averaging the adjacent cell values
!              using the averaging coefficients avgCo.
!
! Input: region = data of current region
!
! Output: fVar = face values of cv.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_fluLesC2F.F90,v 1.5 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FluLesC2F( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: if, ifl, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  INTEGER :: ijkC0, ijkC1, ifg, bcType
  REAL(RFREAL), POINTER :: cv(:,:), avgCo(:,:), fVar(:,:), bfVar(:,:)
  INTEGER, POINTER      :: f2c(:,:), bf2c(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluLesC2F.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FluLesC2F',&
  'TURB_fluLesC2F.F90' )

! get indices and pointers ---------------------------------------------------

  cv         => region%mixt%cv
  f2c        => region%grid%f2c
  avgCo      => region%turb%avgCoI
  fVar       => region%turb%fVar
  bfVar      => region%turb%bfVar
  avgCo(:,:) =  0.5_RFREAL  ! for now, waiting for proper average coefficients

! perform 2-point averaging of rho,rhou1,rhou2 and rhou3 from centers to faces

  DO if = 1,region%grid%nFaces
    ijkC0 = f2c(1,if)
    ijkC1 = f2c(2,if)

    fVar(CV_TURB_DENS,if)=avgCo(2,if)*cv(CV_MIXT_DENS,ijkC0)+ &
                          avgCo(1,if)*cv(CV_MIXT_DENS,ijkC1)
    fVar(CV_TURB_XMOM,if)=avgCo(2,if)*cv(CV_MIXT_XMOM,ijkC0)+ &
                          avgCo(1,if)*cv(CV_MIXT_XMOM,ijkC1)
    fVar(CV_TURB_YMOM,if)=avgCo(2,if)*cv(CV_MIXT_YMOM,ijkC0)+ &
                          avgCo(1,if)*cv(CV_MIXT_YMOM,ijkC1)
    fVar(CV_TURB_ZMOM,if)=avgCo(2,if)*cv(CV_MIXT_ZMOM,ijkC0)+ &
                          avgCo(1,if)*cv(CV_MIXT_ZMOM,ijkC1)
  ENDDO  ! if

  ifg = 0
  DO iPatch = 1,region%grid%nPatches
    patch  => region%patches(iPatch)
    bf2c   => patch%bf2c
    bcType =  patch%bcType

    IF (bcType>=BC_NOSLIPWALL .AND. bcType>=BC_NOSLIPWALL+BC_RANGE) THEN   

      DO ifl = 1,patch%nBFaces
        ijkC0 = bf2c(ifl)
        ifg   = ifg + 1

        bfVar(CV_TURB_DENS,ifg) = cv(CV_MIXT_DENS,ijkC0)
        bfVar(CV_TURB_XMOM,ifg) = 0._RFREAL
        bfVar(CV_TURB_YMOM,ifg) = 0._RFREAL
        bfVar(CV_TURB_ZMOM,ifg) = 0._RFREAL
      ENDDO  ! ifl
    ELSE
      DO ifl = 1,patch%nBFaces
        ijkC0 = bf2c(ifl)
        ifg   = ifg + 1

        bfVar(CV_TURB_DENS,ifg) = cv(CV_MIXT_DENS,ijkC0)
        bfVar(CV_TURB_XMOM,ifg) = cv(CV_MIXT_XMOM,ijkC0)
        bfVar(CV_TURB_YMOM,ifg) = cv(CV_MIXT_YMOM,ijkC0)
        bfVar(CV_TURB_ZMOM,ifg) = cv(CV_MIXT_ZMOM,ijkC0)
      ENDDO  ! ifl
    ENDIF

  ENDDO  ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FluLesC2F

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_fluLesC2F.F90,v $
! Revision 1.5  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/05/28 02:00:01  wasistho
! update unstructured grid LES
!
! Revision 1.2  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.1  2004/03/19 02:54:58  wasistho
! prepared for RFLU
!
!
!******************************************************************************







