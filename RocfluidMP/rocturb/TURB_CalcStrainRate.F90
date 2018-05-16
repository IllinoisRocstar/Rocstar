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
! Purpose: Compute strain rate tensor S_ij from any given velocity gradients.
!
! Description: The strain rate tensors are computed for i, j and k faces.
!
! Input: region  = data of current region
!        ibn,ien = begin and end nodes index
!        grIndx  = ID index of the velocity gradients
!        gradi,j,k  = velocity gradients at i, j and k faces
!
! Output: sRateI,J,K = strain rate tensor at i, j and k faces
!
! Notes: Strain rate values at all nodes are defined.
!        RFLU data structure makes use only of I components arrays.
!
!******************************************************************************
!
! $Id: TURB_CalcStrainRate.F90,v 1.5 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE TURB_CalcStrainRate( region,ibn,ien,grIndx,gradi,gradj,gradk, &
                                sRateI,sRateJ,sRateK )
#endif
#ifdef RFLU
SUBROUTINE TURB_CalcStrainRate( region,ibn,ien,grIndx,gradf,sRateI )
#endif

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ibn, ien
#ifdef RFLO
  INTEGER        :: grIndx(9)
  REAL(RFREAL), POINTER :: gradi(:,:),gradj(:,:),gradk(:,:)
  REAL(RFREAL), POINTER :: sRateI(:,:),sRateJ(:,:),sRateK(:,:)
#endif
#ifdef RFLU
  INTEGER        :: grIndx(3)
  REAL(RFREAL), POINTER :: gradf(:,:,:)
  REAL(RFREAL), POINTER :: sRateI(:,:)
#endif

! ... loop variables
  INTEGER :: iN

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  REAL(RFREAL)          :: oo3,div
  REAL(RFREAL)          :: ux,uy,uz,vx,vy,vz,wx,wy,wz

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_CalcStrainRate.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_CalcStrainRate',&
  'TURB_CalcStrainRate.F90' )

  oo3 = 1.0_RFREAL/3.0_RFREAL

! get strain rate tensor at cell faces --------------------------------------

  DO iN = ibn,ien

! - strain rate tensor at I faces

#ifdef RFLO
    ux  = gradi(grIndx(1),iN)
    vx  = gradi(grIndx(2),iN)
    wx  = gradi(grIndx(3),iN)
    uy  = gradi(grIndx(4),iN)
    vy  = gradi(grIndx(5),iN)
    wy  = gradi(grIndx(6),iN)
    uz  = gradi(grIndx(7),iN)
    vz  = gradi(grIndx(8),iN)
    wz  = gradi(grIndx(9),iN)
#endif
#ifdef RFLU
    ux  = gradf(XCOORD,grIndx(1),iN)
    uy  = gradf(YCOORD,grIndx(1),iN)
    uz  = gradf(ZCOORD,grIndx(1),iN)
    vx  = gradf(XCOORD,grIndx(2),iN)
    vy  = gradf(YCOORD,grIndx(2),iN)
    vz  = gradf(ZCOORD,grIndx(2),iN)
    wx  = gradf(XCOORD,grIndx(3),iN)
    wy  = gradf(YCOORD,grIndx(3),iN)
    wz  = gradf(ZCOORD,grIndx(3),iN)
#endif

    div = oo3*(ux+vy+wz)

    sRateI(E11,iN) = 2.0_RFREAL*(ux-div)
    sRateI(E12,iN) = uy+vx
    sRateI(E13,iN) = uz+wx

    sRateI(E22,iN) = 2.0_RFREAL*(vy-div)
    sRateI(E23,iN) = vz+wy

    sRateI(E33,iN) = 2.0_RFREAL*(wz-div)

#ifdef RFLO
! - strain rate tensor at J faces

    ux  = gradj(grIndx(1),iN)
    vx  = gradj(grIndx(2),iN)
    wx  = gradj(grIndx(3),iN)
    uy  = gradj(grIndx(4),iN)
    vy  = gradj(grIndx(5),iN)
    wy  = gradj(grIndx(6),iN)
    uz  = gradj(grIndx(7),iN)
    vz  = gradj(grIndx(8),iN)
    wz  = gradj(grIndx(9),iN)

    div = oo3*(ux+vy+wz)

    sRateJ(E11,iN) = 2.0_RFREAL*(ux-div)
    sRateJ(E12,iN) = uy+vx
    sRateJ(E13,iN) = uz+wx

    sRateJ(E22,iN) = 2.0_RFREAL*(vy-div)
    sRateJ(E23,iN) = vz+wy

    sRateJ(E33,iN) = 2.0_RFREAL*(wz-div)

! - strain rate tensor at K faces

    ux  = gradk(grIndx(1),iN)
    vx  = gradk(grIndx(2),iN)
    wx  = gradk(grIndx(3),iN)
    uy  = gradk(grIndx(4),iN)
    vy  = gradk(grIndx(5),iN)
    wy  = gradk(grIndx(6),iN)
    uz  = gradk(grIndx(7),iN)
    vz  = gradk(grIndx(8),iN)
    wz  = gradk(grIndx(9),iN)

    div = oo3*(ux+vy+wz)

    sRateK(E11,iN) = 2.0_RFREAL*(ux-div)
    sRateK(E12,iN) = uy+vx
    sRateK(E13,iN) = uz+wx

    sRateK(E22,iN) = 2.0_RFREAL*(vy-div)
    sRateK(E23,iN) = vz+wy

    sRateK(E33,iN) = 2.0_RFREAL*(wz-div)
#endif

  ENDDO     ! iN

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CalcStrainRate

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_CalcStrainRate.F90,v $
! Revision 1.5  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/05/17 20:32:55  wasistho
! reordering gradIndx for more efficient cache memory addressing
!
! Revision 1.2  2004/03/19 02:44:55  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







