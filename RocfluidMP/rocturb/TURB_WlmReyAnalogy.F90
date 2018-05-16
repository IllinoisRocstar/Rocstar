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
! Purpose: Compute wall heat flux (heat transfer) assuming Reynolds analogy 
!          between heat and momentum transfer.
!
! Description: Reynolds analogy: 
!              St/(0.5*Cf) = Prt^(-0.66). 
!              With Stanton number, St = qw/(rho*cp*(Tw-Te)*ue), and
!              assuming (Te-Tw)/ue = (T1-Tw)/u1, heat flux (written dotless) 
!              qw = tau_w*cp*(Tw-T1)/u1. 
!              Heat transfer coefficient: cq = qw/(Tw-Te)
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: wall heat flux.
!
! Notes: - Avoiding complication determining BL-edge, heat transfer coefficient
!          if desired should, currently, be computed externally.
!        - Total wall stress already stored in vals(:,WLM_VALS_HFLUX) to use.
!
!******************************************************************************
!
! $Id: TURB_WlmReyAnalogy.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmReyAnalogy( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: ijBeg, ijEnd, indCp, bcOpt, distrib, ijkVal, ijkC, ijkCe
  REAL(RFREAL)          :: tWall, tauWall, abVel, ratio, prt, cp, utau
  REAL(RFREAL), POINTER :: dv(:,:), gv(:,:), vals(:,:), mVals(:,:)

#ifdef RFLO
  INTEGER :: n1, n2, iOff, iCOff, ijCOff
  INTEGER :: ilev, idir, jdir, kdir, iedge, jedge, kedge, nedge, nEta
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER :: cv(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmReyAnalogy.F90,v $ $Revision: 1.6 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmReyAnalogy',&
  'TURB_WlmReyAnalogy.F90' )

! get dimensions, pointers and parameters -------------------------------------

  bcOpt =  patch%mixt%switches(BCSWI_NOSLIP_ADIABAT)
  mVals => patch%mixt%vals
  vals  => patch%valBola%vals 

#ifdef RFLO
  n1    = ABS(patch%l1end-patch%l1beg)
  n2    = ABS(patch%l2end-patch%l2beg)
  iOff  = n1 + 1
  ijBeg = IndIJ( 0, 0,iOff)
  ijEnd = IndIJ(n1,n2,iOff)
#endif
#ifdef RFLU
  ijBeg = 1
  ijEnd = patch%nBFaces
#endif

! adiabatic wall : skip heat flux computation and set it to zero

  IF (bcOpt == BCOPT_ADIABAT) THEN
    DO ijkVal=ijBeg,ijEnd
      vals(ijkVal,WLM_VALS_HFLUX) = 0._RFREAL
    ENDDO
    GOTO 999
  ENDIF

! isothermal wall: proceed with heat flux computation

#ifdef RFLO
  ilev   =  region%currLevel

  dv     => region%levels(ilev)%mixt%dv
  gv     => region%levels(ilev)%mixt%gv
  indCp  =  region%levels(iLev)%mixt%indCp
  prt    =  region%levels(iLev)%mixt%prTurb
#endif
#ifdef RFLU
  cv     => region%mixt%cv
  dv     => region%mixt%dv
  gv     => region%mixt%gv
  indCp  =  region%mixtInput%indCp
  prt    =  region%mixtInput%prTurb
#endif
  ratio   = prt**(-0.66_RFREAL)
  distrib = patch%mixt%distrib

#ifdef RFLO
  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

  nEta  = ABS(idir)*region%levels(ilev)%grid%ipc + &
          ABS(jdir)*region%levels(ilev)%grid%jpc + &
          ABS(kdir)*region%levels(ilev)%grid%kpc
  nedge = MIN( 8,nEta )
  iedge = idir*nedge
  jedge = jdir*nedge
  kedge = kdir*nedge

! Compute wall heat flux. 
! Note that total wall stress has been stored in vals(:,WLM_VALS_HFLUX).

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkVal = ABS(idir)*IndIJ(j-jbeg,k-kbeg,jend-jbeg+1) + &
                 ABS(jdir)*IndIJ(k-kbeg,i-ibeg,kend-kbeg+1) + &
                 ABS(kdir)*IndIJ(i-ibeg,j-jbeg,iend-ibeg+1)
        ijkC   = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)
        ijkCe  = IndIJK(i+iedge,j+jedge,k+kedge,iCOff,ijCOff)

        abVel  = SQRT( dv(DV_MIXT_UVEL,ijkCe)**2+dv(DV_MIXT_VVEL,ijkCe)**2+ &
                       dv(DV_MIXT_WVEL,ijkCe)**2 )
#endif
#ifdef RFLU
! Specific Rocflu, check the state of cv first
  IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) &
                            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  ibeg  =  1
  iend  =  patch%nBFaces

  DO i=ibeg,iend
        ijkVal = i
        ijkC   = patch%bf2c(i)
        ijkCe  = ijkC

        abVel  = SQRT( cv(CV_MIXT_XMOM,ijkCe)**2+cv(CV_MIXT_YMOM,ijkCe)**2+ &
                       cv(CV_MIXT_ZMOM,ijkCe)**2 )/cv(CV_MIXT_DENS,ijkCe)
#endif
        tWall  = mVals(BCDAT_NOSLIP_TWALL,distrib*ijkVal)
        cp     = gv(GV_MIXT_CP,ijkC*indCp)
        tauWall= vals(ijkVal,WLM_VALS_HFLUX)
        vals(ijkVal,WLM_VALS_HFLUX) = ratio*tauWall*cp* &
                                      (tWall-dv(DV_MIXT_TEMP,ijkCe))/abVel

!wlmCheckprobe-----------------------------------------------------------------
!        utau = SQRT( tauWall/vals(ijkVal,WLM_VALS_DENS) )
!        write(*,*)region%procId,nEta,ijkVal,vals(ijkVal,WLM_VALS_DENS),utau, &
!                  tauWall,vals(ijkVal,WLM_VALS_HFLUX)
!------------------------------------------------------------------------------
#ifdef RFLO
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  ENDDO       ! i
#endif

999 CONTINUE

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmReyAnalogy

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmReyAnalogy.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:40:42  mparmar
! Renamed patch variables
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:01  wasistho
! changed nomenclature
!
! Revision 1.2  2004/03/02 03:50:55  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







