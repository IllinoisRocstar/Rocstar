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
! Purpose: calculate max. allowable local/global time step in the case
!          of viscous flow.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%dt   = local time step
!         regions%levels%srad = convective spectral radii
!         global%dtMin        = global time step for all regions on this
!                               processor(if unsteady flow).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_TimeStepViscous.F90,v 1.3 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_TimeStepViscous( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
        RFLO_GetNodeOffset, RFLO_copyVectorPatches, RFLO_copyVectorEdges, &
        RFLO_copyVectorCorners, RFLO_copyMatrixPatches, RFLO_copyMatrixEdges, &
        RFLO_copyMatrixCorners, RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: indCp, indMol, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend, turbModel
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkN, ijkN1, indSvel

  REAL(RFREAL) :: rgas, cpgas, gamma, f1, f2, fac, fmue, dtv, prLam, prTurb
  REAL(RFREAL) :: rrho, u, v, w, sx, sy, sz, sVel, dS2, dS, vc, cs, sumSrad
  REAL(RFREAL) :: srvi, srvj, srvk, dtMin
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), tv(:,:), dt(:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:), vol(:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)
  REAL(RFREAL), POINTER :: srad(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_TimeStepViscous',&
  'RFLO_TimeStepViscous.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                           kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv => region%levels(iLev)%mixt%cv
  dv => region%levels(iLev)%mixt%dv
  gv => region%levels(iLev)%mixt%gv
  tv => region%levels(iLev)%mixt%tv
  dt => region%levels(iLev)%dt

  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  vol    => region%levels(iLev)%grid%vol
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel

  srad => region%levels(iLev)%mixt%srad

  prLam     = region%levels(iLev)%mixt%prLam
  prTurb    = region%levels(iLev)%mixt%prTurb
  indCp     = region%levels(iLev)%mixt%indCp
  indMol    = region%levels(iLev)%mixt%indMol
  turbModel = region%mixtInput%turbModel

! local time step, spectral radii ---------------------------------------------

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC  = IndIJK(i,j,k,iCoff,ijCOff)
        ijkN  = IndIJK(i,j,k,iNOff,ijNOff)
        rrho  = 1._RFREAL/cv(CV_MIXT_DENS,ijkC)
        u     = cv(CV_MIXT_XMOM,ijkC)*rrho
        v     = cv(CV_MIXT_YMOM,ijkC)*rrho
        w     = cv(CV_MIXT_ZMOM,ijkC)*rrho
        rgas  = 8314.3_RFREAL/gv(GV_MIXT_MOL,ijkC*indMol)
        cpgas = gv(GV_MIXT_CP,ijkC*indCp)
        gamma = cpgas/(cpgas-rgas)
        f1    = 4._RFREAL*rrho/3._RFREAL
        f2    = gamma*rrho
        fac   = MAX(f1,f2)
        IF (turbModel == TURB_MODEL_NONE) THEN
          fmue = tv(TV_MIXT_MUEL,ijkC)/prLam
        ELSE
          fmue = tv(TV_MIXT_MUEL,ijkC)/prLam + tv(TV_MIXT_MUET,ijkC)/prTurb
        ENDIF
        dtv   = fac*fmue/vol(ijkC)

        ijkN1 = IndIJK(i+1,j,k,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,ijkN1))
        sVel  = 0.5_RFREAL*(siVel(ijkN*indSvel)+siVel(ijkN1*indSvel))
        dS2   = sx*sx + sy*sy + sz*sz
        dS    = SQRT(dS2)
        vc    = sx*u + sy*v + sz*w - sVel
        cs    = dv(DV_MIXT_SOUN,ijkC)*dS
        srvi  = dtv*dS2
        srad(ICOORD,ijkC) = ABS(vc) + cs

        ijkN1 = IndIJK(i,j+1,k,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijkN1))
        sVel  = 0.5_RFREAL*(sjVel(ijkN*indSvel)+sjVel(ijkN1*indSvel))
        dS2   = sx*sx + sy*sy + sz*sz
        dS    = SQRT(dS2)
        vc    = sx*u + sy*v + sz*w - sVel
        cs    = dv(DV_MIXT_SOUN,ijkC)*dS
        srvj  = dtv*dS2
        srad(JCOORD,ijkC) = ABS(vc) + cs

        ijkN1 = IndIJK(i,j,k+1,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(sk(XCOORD,ijkN)+sk(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(sk(YCOORD,ijkN)+sk(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(sk(ZCOORD,ijkN)+sk(ZCOORD,ijkN1))
        sVel  = 0.5_RFREAL*(skVel(ijkN*indSvel)+skVel(ijkN1*indSvel))
        dS2   = sx*sx + sy*sy + sz*sz
        dS    = SQRT(dS2)
        vc    = sx*u + sy*v + sz*w - sVel
        cs    = dv(DV_MIXT_SOUN,ijkC)*dS
        srvk  = dtv*dS2
        srad(KCOORD,ijkC) = ABS(vc) + cs

        sumSrad  = srad(ICOORD,ijkC) + srad(JCOORD,ijkC) + srad(KCOORD,ijkC) &
                   + 4._RFREAL*(srvi+srvj+srvk)
        sumSrad  = MAX(sumSrad,1.E-30_RFREAL)
        dt(ijkC) = vol(ijkC)/sumSrad
      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! treat dummy cells -----------------------------------------------------------

  CALL RFLO_copyVectorPatches( iLev,region,dt )
  CALL RFLO_copyVectorEdges( iLev,region,dt )
  CALL RFLO_copyVectorCorners( iLev,region,dt )

  CALL RFLO_copyMatrixPatches( iLev,region,srad )
  CALL RFLO_copyMatrixEdges( iLev,region,srad )
  CALL RFLO_copyMatrixCorners( iLev,region,srad )

! global time step (this processor) -------------------------------------------

  IF (global%flowType == FLOW_UNSTEADY) THEN
    dtMin = global%dtMin

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          ijkC  = IndIJK(i,j,k,iCoff,ijCOff)
          dtMin = MIN(dtMin,dt(ijkC))
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

    global%dtMin = dtMin
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_TimeStepViscous

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_TimeStepViscous.F90,v $
! Revision 1.3  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.17  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.12  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.11  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/09/02 23:09:47  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.9  2002/08/30 18:25:55  jblazek
! Forgot to multiply grid speeds by face area ...
!
! Revision 1.8  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.7  2002/08/29 22:45:56  jblazek
! Added support for moving grids.
!
! Revision 1.6  2002/08/28 00:03:56  wasistho
! mueT put back to consist only turbulent part
!
! Revision 1.5  2002/08/24 03:17:15  wasistho
! put safety within #ifdef TURB
!
! Revision 1.4  2002/08/23 03:18:37  wasistho
! due to mueT as Total instead of Turbulent
!
! Revision 1.3  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
!******************************************************************************







