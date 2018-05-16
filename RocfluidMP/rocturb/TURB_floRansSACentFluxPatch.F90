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
! Purpose: compute central convective fluxe of SA eq. through a patch
!          by using an average of variables.
!
! Description: none.
!
! Input: region = data of current region
!        patch  = current patch.
!
! Output: region%levels%turb%rhs = convective fluxes added to SA residual.
!
! Notes: grid speeds already contain face area.
!
!******************************************************************************
!
! $Id: TURB_floRansSACentFluxPatch.F90,v 1.7 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloRansSACentFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetNodeOffset, BcondInjectionPerf
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, n1, n2

! ... local variables
  INTEGER :: iLev, lbound, bcType, distrib, gasModel, indCp, indMol
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, i2d, nOff
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCD, ijkCB0, ijkCB1, ijkNB
  INTEGER :: inode, jnode, knode, indSvel, acId0, acId1

  REAL(RFREAL)          :: sgn, rhoa, rhoua, rhova, rhowa, rhoea, rhoVrel(3)
  REAL(RFREAL)          :: rnutila, pa, mRate, tBurn, dS, sf(3)
  REAL(RFREAL)          :: ua, va, wa, uinj, vinj, winj, vcont, sv, ac0, ac1
  REAL(RFREAL), POINTER :: dv(:,:), gv(:,:), tcv(:,:), trhs(:,:)
  REAL(RFREAL), POINTER :: avgCo(:,:), sFace(:,:), sVel(:), vals(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_FloRansSACentFluxPatch',&
  'TURB_floRansSACentFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  bcType    = patch%bcType
  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  gasModel = region%mixtInput%gasModel
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  indSvel   = region%levels(iLev)%grid%indSvel

  dv   => region%levels(iLev)%mixt%dv
  gv   => region%levels(iLev)%mixt%gv
  tcv  => region%levels(iLev)%turb%cv
  trhs => region%levels(iLev)%turb%rhs
  vals => patch%mixt%vals

! to take the right face vector and make it point outwards

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

! get the appropriate face vector and grid speed

  IF (lbound==1 .OR. lbound==2) THEN
    avgCo => region%levels(iLev)%grid%c2fCoI
    sFace => region%levels(iLev)%grid%si
    sVel  => region%levels(iLev)%grid%siVel
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(iLev)%grid%sj
    sVel  => region%levels(iLev)%grid%sjVel
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(iLev)%grid%sk
    sVel  => region%levels(iLev)%grid%skVel
  ENDIF

! stationary grid -------------------------------------------------------------

  IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) THEN

  ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

  ELSE IF (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE) THEN

  ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

! injection boundary (as wall if mass flow rate <= 0)

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkCB1 = IndIJK(i+idir ,j+jdir ,k+kdir ,iCOff,ijCOff)  ! interior
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          ac0    = avgCo(acId0,ijkNB)
          ac1    = avgCo(acId1,ijkNB)
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

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
          i2d        = distrib * IndIJ(n1,n2,nOff)
          mRate      = vals(BCDAT_INJECT_MFRATE,i2d)
          tBurn      = vals(BCDAT_INJECT_TEMP  ,i2d)
          rhoVrel(1) = vals(BCDAT_INJECT_RFVFU ,i2d)
          rhoVrel(2) = vals(BCDAT_INJECT_RFVFV ,i2d)
          rhoVrel(3) = vals(BCDAT_INJECT_RFVFW ,i2d)

          IF (mRate > 0._RFREAL) THEN        ! surface burning
!            rnutila = ac0*tcv(CV_SA_NUTIL,ijkCB0)+ &
!                      ac1*tcv(CV_SA_NUTIL,ijkCD)
            rnutila = 0.5_RFREAL*(3._RFREAL*tcv(CV_SA_NUTIL,ijkCB0)- &
                                            tcv(CV_SA_NUTIL,ijkCB1))
            dS     = SQRT(sf(1)*sf(1)+sf(2)*sf(2)+sf(3)*sf(3))

            IF (gasModel == GAS_MODEL_TCPERF) THEN
              CALL BcondInjectionPerf( distrib,mRate,tBurn,rhoVrel, &
                                       sf(1)/dS,sf(2)/dS,sf(3)/dS, &
                                       gv(GV_MIXT_CP  ,ijkCB0*indCp ), &
                                       gv(GV_MIXT_MOL ,ijkCB0*indMol), &
                                       dv(DV_MIXT_PRES,ijkCB0       ), &
                                       rhoa,rhoua,rhova,rhowa,rhoea,pa, &
                                       uinj,vinj,winj )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,__LINE__ )
            ENDIF
            vcont = uinj*sf(1) + vinj*sf(2) + winj*sf(3)

          ELSE                               ! not burning - slip/noslip wall
            rnutila = 0._RFREAL
            vcont   = 0._RFREAL
          ENDIF

          trhs(CV_SA_NUTIL,ijkCB0) = trhs(CV_SA_NUTIL,ijkCB0) + &
                                     vcont*rnutila

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! non-conforming region interface

  ELSE IF (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! everything else

  ELSE

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
          ac0    = avgCo(acId0,ijkNB)
          ac1    = avgCo(acId1,ijkNB)

          rnutila= ac0*tcv(CV_SA_NUTIL,ijkCB0)+ac1*tcv(CV_SA_NUTIL,ijkCD)
          ua     = ac0*dv(DV_MIXT_UVEL,ijkCB0)+ac1*dv(DV_MIXT_UVEL,ijkCD)
          va     = ac0*dv(DV_MIXT_VVEL,ijkCB0)+ac1*dv(DV_MIXT_VVEL,ijkCD)
          wa     = ac0*dv(DV_MIXT_WVEL,ijkCB0)+ac1*dv(DV_MIXT_WVEL,ijkCD)
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
          sv     = sVel(ijkNB*indSvel)
          vcont  = ua*sf(1)+va*sf(2)+wa*sf(3) - sv

          trhs(CV_SA_NUTIL,ijkCB0) = trhs(CV_SA_NUTIL,ijkCB0) + &
                                     vcont*rnutila
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  ENDIF         ! bcType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_FloRansSACentFluxPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floRansSACentFluxPatch.F90,v $
! Revision 1.7  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:40:49  mparmar
! Renamed patch variables
!
! Revision 1.4  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.3  2004/08/02 21:55:11  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.3  2003/10/21 22:08:47  wasistho
! modified grid speed as it already contain face area
!
! Revision 1.2  2003/10/20 20:27:40  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







