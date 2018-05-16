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
!          through a patch
!
! Description: this routine works in the same way as ViscousFluxEddy
!              but applied on region patches
!
! Input: region  = data of current region.
!        patch   = current patch.
!        indxMu  = TV component index for dynamic viscosity Mu 
!        indxTCo = TV component index for thermal conductivity TCo 
!        tv      = transport variables, mu and kappa
!
! Output: region%levels%mixt%diss = viscous fluxes added to the dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ViscousFluxPatch.F90,v 1.4 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ViscousFluxPatch( region,patch,indxMu,indxTCo,tv )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"

  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch
  INTEGER        :: indxMu, indxTCo
  REAL(RFREAL), POINTER :: tv(:,:)

! ... loop variables
  INTEGER :: i, j, k, iC

! ... local variables
  INTEGER :: ilev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inode, jnode, knode, idir, jdir, kdir
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCB0, ijkCD, ijkNB, acId0, acId1
  INTEGER :: n1, n2, nOff, j2d, aeroCoeff

  REAL(RFREAL)          :: oo3, beta, muf, kpaf, velf(3), div
  REAL(RFREAL)          :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty, tz, &
                           tgradf, fd(4)
  REAL(RFREAL)          :: sgn, sf(3), sij(3,3), ac0, ac1
  REAL(RFREAL)          :: rRef, vRef, rCfRef, rChRef
  REAL(RFREAL), POINTER :: avgCo(:,:), diss(:,:)

  REAL(RFREAL), POINTER :: dv(:,:), grad(:,:), sFace(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ViscousFluxPatch',&
  'RFLO_ViscousFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  bcType = patch%bcType

  ilev   = region%currLevel
  lbound = patch%lbound
  
  dv   => region%levels(ilev)%mixt%dv
  diss => region%levels(ilev)%mixt%diss  

! get coefficients -----------------------------------------------------------

  oo3  = 1.0_RFREAL/3.0_RFREAL
  beta = region%mixtInput%betrk(region%irkStep)

  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  aeroCoeff = region%global%aeroCoeffs
  rRef      = region%global%refDensity
  vRef      = region%global%refVelocity
  
  rCfRef = 2.0_RFREAL/(rRef*vRef*vRef)
  rChRef = 2.0_RFREAL/(rRef*vRef*vRef*vRef)

! take the right face vector and make it point outwards ----------------------

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  acId0 = 2
  acId1 = 1
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
    acId0 = 1
    acId1 = 2
  ENDIF

! get the appropriate face vector --------------------------------------------

  IF (lbound==1 .OR. lbound==2) THEN
    avgCo => region%levels(iLev)%grid%c2fCoI
    sFace => region%levels(ilev)%grid%si
    grad  => region%levels(ilev)%mixt%gradi
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(ilev)%grid%sj
    grad  => region%levels(ilev)%mixt%gradj
  ELSE
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(ilev)%grid%sk
    grad  => region%levels(ilev)%mixt%gradk
  ENDIF

! non-conforming region interface --------------------------------------------

  IF (bcType>=BC_regionINT .AND. bcType<=BC_regionINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_regNONCONF .AND. bcType<=BC_regNONCONF+BC_RANGE) THEN

! everything else

  ELSE

! flux in the direction normal to the patch ----------------------------------

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! bnd cells
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)  ! bnd nodes
          ac0    = avgCo(acId0,ijkNB)
          ac1    = avgCo(acId1,ijkNB)

          muf  = ac0*tv(TV_MIXT_MUEL,ijkCB0)+ac1*tv(TV_MIXT_MUEL,ijkCD)
          kpaf = ac0*tv(TV_MIXT_TCOL,ijkCB0)+ac1*tv(TV_MIXT_TCOL,ijkCD)

          velf(1)= ac0*dv(DV_MIXT_UVEL,ijkCB0)+ac1*dv(DV_MIXT_UVEL,ijkCD)
          velf(2)= ac0*dv(DV_MIXT_VVEL,ijkCB0)+ac1*dv(DV_MIXT_VVEL,ijkCD)
          velf(3)= ac0*dv(DV_MIXT_WVEL,ijkCB0)+ac1*dv(DV_MIXT_WVEL,ijkCD)
                             
          ux  = grad(GR_MIXT_UX,ijkNB)
          uy  = grad(GR_MIXT_UY,ijkNB)
          uz  = grad(GR_MIXT_UZ,ijkNB)
          vx  = grad(GR_MIXT_VX,ijkNB)
          vy  = grad(GR_MIXT_VY,ijkNB)
          vz  = grad(GR_MIXT_VZ,ijkNB)
          wx  = grad(GR_MIXT_WX,ijkNB)
          wy  = grad(GR_MIXT_WY,ijkNB)
          wz  = grad(GR_MIXT_WZ,ijkNB)
          tx  = grad(GR_MIXT_TX,ijkNB)
          ty  = grad(GR_MIXT_TY,ijkNB)
          tz  = grad(GR_MIXT_TZ,ijkNB)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

          tgradf= (tx*sf(1)+ty*sf(2)+tz*sf(3))

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

          fd(1) = muf*(sij(1,1)*sf(1)+sij(1,2)*sf(2)+sij(1,3)*sf(3))
          fd(2) = muf*(sij(2,1)*sf(1)+sij(2,2)*sf(2)+sij(2,3)*sf(3))
          fd(3) = muf*(sij(3,1)*sf(1)+sij(3,2)*sf(2)+sij(3,3)*sf(3))
          fd(4) = DOT_PRODUCT(fd(1:3),velf(1:3)) + kpaf*tgradf

          diss(CV_MIXT_XMOM,ijkCB0) = diss(CV_MIXT_XMOM,ijkCB0)+fd(1)*beta
          diss(CV_MIXT_YMOM,ijkCB0) = diss(CV_MIXT_YMOM,ijkCB0)+fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkCB0) = diss(CV_MIXT_ZMOM,ijkCB0)+fd(3)*beta
          diss(CV_MIXT_ENER,ijkCB0) = diss(CV_MIXT_ENER,ijkCB0)+fd(4)*beta

! ------- Set friction and heat-transfer coefficients

          IF      (lbound==1 .OR. lbound==2) THEN
            n1 = j - jbeg
            n2 = k - kbeg
          ELSE IF (lbound==3 .OR. lbound==4) THEN
            n1 = k - kbeg
            n2 = i - ibeg
          ELSE IF (lbound==5 .OR. lbound==6) THEN
            n1 = i - ibeg
            n2 = j - jbeg
          ENDIF
          j2d  = aeroCoeff * IndIJ(n1,n2,nOff)
          patch%cf(XCOORD,j2d) = rCfRef*fd(1)
          patch%cf(YCOORD,j2d) = rCfRef*fd(2)
          patch%cf(ZCOORD,j2d) = rCfRef*fd(3)
          patch%ch(       j2d) = rChRef*fd(4)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
  ENDIF 

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ViscousFluxPatch

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_ViscousFluxPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/13 03:42:28  wasistho
! defined aero coeffs
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/02 23:13:59  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.15  2003/12/04 03:23:09  haselbac
! Removed RFLU code segments, moved to RFLU_ViscousFluxesPatches
!
! Revision 1.14  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.10  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.9  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.8  2002/09/09 14:08:12  haselbac
! mixtInput now under regions, added some things, corrected a few bugs
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/08/01 01:40:03  wasistho
! Made generic
!
! Revision 1.4  2002/07/27 08:13:52  wasistho
! prepared for rocturb implementation
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:38:48  wasistho
! made compliant with CODING RULE
!
!******************************************************************************







