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
! Purpose: compute viscous fluxes based on viscosity principle: mu*S_ij
!
! Description: this routine compute S_ij, while mu is given
!              mu can be mu_l or mu_tot=mu_l+mu_t with l=laminar, t=turbulent
!
! Input: region  = data of current region
!        indxMu  = TV component index for dynamic viscosity Mu 
!        indxTCo = TV component index for thermal conductivity TCo 
!        tv      = transport variables, mu and kappa
!
! Output: region%levels%mixt%diss = viscous fluxes added to the dissipation.
!
! Notes: indxMu can have value TV_MIXT_MUEL (laminar) or TV_MIXT_MUET (total)
!       indxTCo can have value TV_MIXT_TCOL (laminar) or TV_MIXT_TCOT (total)
!
!******************************************************************************
!
! $Id: RFLO_ViscousFlux.F90,v 1.3 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ViscousFlux( region,indxMu,indxTCo,tv )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE ModInterfaces, ONLY : RFLO_ViscousFluxPatch

#include "Indexing.h"
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: indxMu, indxTCo
  REAL(RFREAL), POINTER :: tv(:,:)

! ... loop variables
  INTEGER :: i, j, k, iC, ipatch

! ... local variables
  INTEGER            :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER            :: ilev, iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkN
  INTEGER, PARAMETER :: IDIR=1, JDIR=2, KDIR=3

  REAL(RFREAL)          :: oo3, beta, muf, kpaf, velf(3), div
  REAL(RFREAL)          :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz, &
                           tgradf, fd(4)
  REAL(RFREAL)          :: sFace(3), sij(3,3)
  REAL(RFREAL), POINTER :: avgCo(:,:), sf(:,:), diss(:,:)
  REAL(RFREAL), POINTER :: dv(:,:), grad(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ViscousFlux',&
  'RFLO_ViscousFlux.F90' )

! get dimensions and pointers ------------------------------------------------

  ilev =  region%currLevel
  dv   => region%levels(ilev)%mixt%dv
  diss => region%levels(ilev)%mixt%diss
  
  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

! get coefficients -----------------------------------------------------------

  oo3  = 1.0_RFREAL/3.0_RFREAL
  beta = region%mixtInput%betrk(region%irkStep)

! interior fluxes -------------------------------------------------------------

  CALL ComputeFlux( IDIR )
  CALL ComputeFlux( JDIR )
  CALL ComputeFlux( KDIR )

! fluxes through boundaries ---------------------------------------------------

  DO ipatch=1,region%nPatches
    CALL RFLO_ViscousFluxPatch( region,region%levels(ilev)%patches(ipatch), &
                                indxMu,indxTCo,tv )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! =============================================================================
!   Flux computation subroutines
! =============================================================================

CONTAINS

  SUBROUTINE ComputeFlux( ijk )

! ... parameters
    INTEGER   :: ijk

! ... local variables
    INTEGER   :: ibeg,iend,jbeg,jend,kbeg,kend, iadd,jadd,kadd
    REAL(RFREAL) :: ac0, ac1

! - Set limits and pointers ---------------------------------------------------

    IF (ijk==IDIR) THEN
      ibeg = ipcbeg+1
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = -1
      jadd = 0
      kadd = 0
      grad  => region%levels(ilev)%mixt%gradi
      sf    => region%levels(ilev)%grid%si
      avgCo => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==JDIR) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg+1
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = 0
      jadd = -1
      kadd = 0
      grad  => region%levels(ilev)%mixt%gradj
      sf    => region%levels(ilev)%grid%sj
      avgCo => region%levels(iLev)%grid%c2fCoJ
    ELSEIF (ijk==KDIR) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg+1
      kend = kpcend
      iadd = 0
      jadd = 0
      kadd = -1
      grad  => region%levels(ilev)%mixt%gradk
      sf    => region%levels(ilev)%grid%sk
      avgCo => region%levels(iLev)%grid%c2fCoK
    ENDIF

! -- flux in ijk-direction (except through boundary) --------------------------

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkC1 = IndIJK(i+iadd,j+jadd,k+kadd,iCOff,ijCOff)
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ac0   = avgCo(2,ijkN)
          ac1   = avgCo(1,ijkN)
          sFace(1)= sf(XCOORD,ijkN)
          sFace(2)= sf(YCOORD,ijkN)
          sFace(3)= sf(ZCOORD,ijkN)

          muf    = ac0*tv(indxMu      ,ijkC0)+ac1*tv(indxMu      ,ijkC1)
          kpaf   = ac0*tv(indxTCo     ,ijkC0)+ac1*tv(indxTCo     ,ijkC1)
                    
          velf(1)= ac0*dv(DV_MIXT_UVEL,ijkC0)+ac1*dv(DV_MIXT_UVEL,ijkC1)
          velf(2)= ac0*dv(DV_MIXT_VVEL,ijkC0)+ac1*dv(DV_MIXT_VVEL,ijkC1)
          velf(3)= ac0*dv(DV_MIXT_WVEL,ijkC0)+ac1*dv(DV_MIXT_WVEL,ijkC1)
          ux  = grad(GR_MIXT_UX,ijkN)
          uy  = grad(GR_MIXT_UY,ijkN)
          uz  = grad(GR_MIXT_UZ,ijkN)
          vx  = grad(GR_MIXT_VX,ijkN)
          vy  = grad(GR_MIXT_VY,ijkN)
          vz  = grad(GR_MIXT_VZ,ijkN)
          wx  = grad(GR_MIXT_WX,ijkN)
          wy  = grad(GR_MIXT_WY,ijkN)
          wz  = grad(GR_MIXT_WZ,ijkN)
          tx  = grad(GR_MIXT_TX,ijkN)
          ty  = grad(GR_MIXT_TY,ijkN)
          tz  = grad(GR_MIXT_TZ,ijkN)

          tgradf= (tx*sFace(1)+ty*sFace(2)+tz*sFace(3))

          div = oo3*(ux+vy+wz)

          sij(1,1) = 2.0_RFREAL*(ux-div)
          sij(1,2) = uy+vx
          sij(1,3) = uz+wx

          sij(2,1) = sij(1,2)
          sij(2,2) = 2.0_RFREAL*(vy-div)
          sij(2,3) = vz+wy

          sij(3,1) = sij(1,3)
          sij(3,2) = sij(2,3)
          sij(3,3) = 2.0_RFREAL*(wz-div)

          fd(1) = muf*(sij(1,1)*sFace(1)+sij(1,2)*sFace(2)+sij(1,3)*sFace(3))
          fd(2) = muf*(sij(2,1)*sFace(1)+sij(2,2)*sFace(2)+sij(2,3)*sFace(3))
          fd(3) = muf*(sij(3,1)*sFace(1)+sij(3,2)*sFace(2)+sij(3,3)*sFace(3))
          fd(4) = DOT_PRODUCT(fd(1:3),velf(1:3)) + kpaf*tgradf

          diss(CV_MIXT_XMOM,ijkC0) = diss(CV_MIXT_XMOM,ijkC0) + fd(1)*beta
          diss(CV_MIXT_YMOM,ijkC0) = diss(CV_MIXT_YMOM,ijkC0) + fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkC0) = diss(CV_MIXT_ZMOM,ijkC0) + fd(3)*beta
          diss(CV_MIXT_ENER,ijkC0) = diss(CV_MIXT_ENER,ijkC0) + fd(4)*beta

          diss(CV_MIXT_XMOM,ijkC1) = diss(CV_MIXT_XMOM,ijkC1) - fd(1)*beta
          diss(CV_MIXT_YMOM,ijkC1) = diss(CV_MIXT_YMOM,ijkC1) - fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkC1) = diss(CV_MIXT_ZMOM,ijkC1) - fd(3)*beta
          diss(CV_MIXT_ENER,ijkC1) = diss(CV_MIXT_ENER,ijkC1) - fd(4)*beta
          
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE ComputeFlux

END SUBROUTINE RFLO_ViscousFlux

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_ViscousFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/02 23:13:59  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.14  2003/12/04 03:23:07  haselbac
! Removed RFLU code segments, moved to RFLU_ViscousFluxes
!
! Revision 1.13  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.10  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.9  2002/10/27 18:49:21  haselbac
! Removed tabs
!
! Revision 1.8  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.7  2002/09/09 14:07:23  haselbac
! mixtInput now under regions, added patch routine, bug fix for sf
!
! Revision 1.6  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/08/01 01:39:24  wasistho
! Included RFLU capability and made generic
!
! Revision 1.4  2002/07/27 08:13:37  wasistho
! prepared for rocturb implementation
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:38:42  wasistho
! made compliant with CODING RULE
!
!******************************************************************************







