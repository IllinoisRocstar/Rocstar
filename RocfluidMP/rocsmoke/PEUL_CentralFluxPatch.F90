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
! Purpose: compute central convective fluxes for smoke through a patch
!          by using an average of variables.
!
! Description: none.
!
! Input: region = data of current region
!        patch  = current patch.
!
! Output: region%levels%peul%rhs = convective fluxes added to the residual.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_CentralFluxPatch.F90,v 1.5 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_CentralFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  USE INRT_ModParameters

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_patch),  INTENT(IN)    :: patch

! ... loop variables
  INTEGER :: i,j,k,ipt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev,lobound,bcType,nPtypes,n1,n2
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend,idir,jdir,kdir,nOff
  INTEGER :: iCOff,ijCOff,iNOff,ijNOff,ijkCD,ijkCB0,ijkNB
  INTEGER :: inode,jnode,knode,indCp,indMol
  INTEGER :: mixtDistrib,peulDistrib,mixtI2d,peulI2d,gasModel

  REAL(RFREAL) :: gRhoa,gRhoua,gRhova,gRhowa,gRhoea
  REAL(RFREAL) :: sgn,sf(3),dS,vcont,mRate
  REAL(RFREAL) :: tBurn,rhoVrel(3),pa,uinj,vinj,winj,vf
  REAL(RFREAL), POINTER :: sCv(:,:),gCv(:,:),gDv(:,:),gGv(:,:)
  REAL(RFREAL), POINTER :: sRhs(:,:)
  REAL(RFREAL), POINTER :: sFace(:,:),mixtVals(:,:),peulVals(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_CentralFluxPatch.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_CentralFluxPatch',&
  'PEUL_CentralFluxPatch.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers -------------------------------------------------

  iLev    = region%currLevel
  lobound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  bcType      = patch%bcType
  nOff        = ABS(patch%l1end-patch%l1beg) + 1
  mixtDistrib = patch%mixt%distrib
  peulDistrib = patch%peul%distrib
  indCp       = region%levels(iLev)%mixt%indCp
  indMol      = region%levels(iLev)%mixt%indMol
  gasModel   = region%mixtInput%gasModel

  gCv  => region%levels(iLev)%mixt%cv
  gDv  => region%levels(iLev)%mixt%dv
  gGv  => region%levels(iLev)%mixt%gv
  sCv  => region%levels(iLev)%peul%cv
  sRhs => region%levels(iLev)%peul%rhs

  nPtypes = region%peulInput%nPtypes
  IF (region%levels(iLev)%peul%nCv /= nPtypes) &
    CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )

  mixtVals => patch%mixt%vals

  NULLIFY(peulVals)
  IF (ASSOCIATED(patch%peul%vals)) THEN
    peulVals => patch%peul%vals
    IF (LBOUND(peulVals,1) /= 1 .OR. UBOUND(peulVals,1) /= nPtypes) &
      CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )
  ENDIF

! to take the right face vector and make it point outwards

  sgn   = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lobound==2 .OR. lobound==4 .OR. lobound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector and grid speed

  SELECT CASE (lobound)
  CASE (1:2)
    sFace => region%levels(iLev)%grid%si
  CASE (3:4)
    sFace => region%levels(iLev)%grid%sj
  CASE (5:6)
    sFace => region%levels(iLev)%grid%sk
  END SELECT

! stationary grid -------------------------------------------------------------

  SELECT CASE(bcType)

! no mass flux for these cases

  CASE (BC_SLIPWALL  :BC_SLIPWALL  +BC_RANGE, &
        BC_NOSLIPWALL:BC_NOSLIPWALL+BC_RANGE, &
        BC_SYMMETRY  :BC_SYMMETRY  +BC_RANGE)

! injection boundary (no mass flux if mass flow rate <= 0)

  CASE (BC_INJECTION:BC_INJECTION+BC_RANGE)

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          SELECT CASE (lobound)

          CASE (1:2)
            n1 = j - jbeg
            n2 = k - kbeg
          CASE (3:4)
            n1 = k - kbeg
            n2 = i - ibeg
          CASE (5:6)
            n1 = i - ibeg
            n2 = j - jbeg

          END SELECT ! lobound

          mixtI2d = mixtDistrib * IndIJ(n1,n2,nOff)
          mRate   = mixtVals(BCDAT_INJECT_MFRATE,mixtI2d)

! ------- only need to augment source terms when there is a mass flux
! ------- (in contrast to gas case, where there is always a pressure term)

          IF (mRate > 0._RFREAL) THEN

            peulI2d = peulDistrib * IndIJ(n1,n2,nOff)
            ijkCB0  = IndIJK(i,j,k,iCOff,ijCOff)
            ijkCD   = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)
            ijkNB   = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

            sf(1) = sgn*sFace(XCOORD,ijkNB)
            sf(2) = sgn*sFace(YCOORD,ijkNB)
            sf(3) = sgn*sFace(ZCOORD,ijkNB)

            tBurn      = mixtVals(BCDAT_INJECT_TEMP  ,mixtI2d)
            rhoVrel(1) = mixtVals(BCDAT_INJECT_RFVFU ,mixtI2d)
            rhoVrel(2) = mixtVals(BCDAT_INJECT_RFVFV ,mixtI2d)
            rhoVrel(3) = mixtVals(BCDAT_INJECT_RFVFW ,mixtI2d)

            dS = SQRT(sf(1)*sf(1)+sf(2)*sf(2)+sf(3)*sf(3))
            IF (gasModel == GAS_MODEL_TCPERF) THEN
              CALL BcondInjectionPerf( mixtDistrib,mRate,tBurn,rhoVrel,   &
                                       sf(1)/dS,sf(2)/dS,sf(3)/dS,        &
                                       gGv(GV_MIXT_CP  ,ijkCB0*indCp ),   &
                                       gGv(GV_MIXT_MOL ,ijkCB0*indMol),   &
                                       gDv(DV_MIXT_PRES,ijkCB0       ),   &
                                       gRhoa,gRhoua,gRhova,gRhowa,gRhoea, &
                                       pa,uinj,vinj,winj )
            ELSE
              CALL ErrorStop( region%global,ERR_UNKNOWN_BC,__LINE__ )
            ENDIF ! gasModel

            vcont = uinj*sf(1) + vinj*sf(2) + winj*sf(3)

            DO ipt=1,nPtypes
              vf = vcont*peulVals(ipt,peulI2d) ! i.e., vf = vcont*sRhoa/gRhoa
              sRhs(ipt,ijkCB0) = sRhs(ipt,ijkCB0) + vf*gRhoa
            END DO ! ipt

          END IF ! mRate
        ENDDO    ! i
      ENDDO      ! j
    ENDDO        ! k

! do not know how to handle these yet

  CASE (BC_REGIONINT :BC_REGIONINT +BC_RANGE)

  CASE (BC_REGNONCONF:BC_REGNONCONF+BC_RANGE)

! everything else

  CASE DEFAULT

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

          gRhoa  = 0.5_RFREAL*(gCv(CV_MIXT_DENS,ijkCB0) + &
                               gCv(CV_MIXT_DENS,ijkCD))
          gRhoua = 0.5_RFREAL*(gCv(CV_MIXT_XMOM,ijkCB0) + &
                               gCv(CV_MIXT_XMOM,ijkCD))
          gRhova = 0.5_RFREAL*(gCv(CV_MIXT_YMOM,ijkCB0) + &
                               gCv(CV_MIXT_YMOM,ijkCD))
          gRhowa = 0.5_RFREAL*(gCv(CV_MIXT_ZMOM,ijkCB0) + &
                               gCv(CV_MIXT_ZMOM,ijkCD))

          sf(1) = sgn*sFace(XCOORD,ijkNB)
          sf(2) = sgn*sFace(YCOORD,ijkNB)
          sf(3) = sgn*sFace(ZCOORD,ijkNB)

          vcont = (gRhoua*sf(1)+gRhova*sf(2)+gRhowa*sf(3))/gRhoa

          DO ipt=1,nPtypes
            sRhs(ipt,ijkCB0) = sRhs(ipt,ijkCB0) + &
              0.5_RFREAL*vcont*(sCv(ipt,ijkCB0)+sCv(ipt,ijkCD))
          END DO ! ipt

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SELECT    ! bcType

! moving grid (add grid speeds) -----------------------------------------------

  IF (region%mixtInput%moveGrid) THEN
    CALL ErrorStop ( global,ERR_PEUL_MOVEGRID,__LINE__ )
  ENDIF         ! moveGrid

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_CentralFluxPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_CentralFluxPatch.F90,v $
! Revision 1.5  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:40:20  mparmar
! Renamed patch variables
!
! Revision 1.2  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 21:09:27  haselbac
! Initial revision after changing case
!
! Revision 1.8  2004/07/28 15:42:13  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2003/09/25 15:42:57  jferry
! Added mixture source terms due to active smoke (for Detangle interaction)
!
! Revision 1.5  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.4  2003/05/05 17:26:57  jferry
! fixed bug in setting n1, n2
!
! Revision 1.3  2003/04/09 15:11:31  jferry
! removed dependence on temporary storage structures
!
! Revision 1.2  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







