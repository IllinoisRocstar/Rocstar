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
! Purpose: compute SA viscous flux: nutot*d_j(tilde[nu])
!
! Description: this routine compute d_j(tilde[nu]), while nutot is given
!              nutot = nu_l+nu_t with l=laminar, t=turbulent
!
! Input: region  = data of current region
!
! Output: region%levels%turb%diss = viscous flux added to SA dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansSAVisFlux.F90,v 1.12 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansSAVisFlux( region )

  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, RFLO_CalcGradVector

#include "Indexing.h"
#endif
#ifdef RFLU
  USE RFLU_ModDifferentiationFaces, ONLY: RFLU_ComputeGradFacesWrapper
  USE RFLU_ModDifferentiationBFaces, ONLY: RFLU_ComputeGradBFacesWrapper
  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace
#endif
  USE TURB_ModInterfaces, ONLY : TURB_RansSAVisFluxPatch
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region) :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif

! ... loop variables
  INTEGER :: i, j, k, iC, ipatch

! ... local variables
  INTEGER :: iBegV, iEndV, iBegG, iEndG, ijkN, ijkC0, ijkC1
  REAL(RFREAL) :: cb2, opcb2, rSigma, beta, rhoa, mua, rnua
  REAL(RFREAL) :: nuf, nuc, nutilX, nutilY, nutilZ, fd
  REAL(RFREAL) :: sFace(3)
  REAL(RFREAL), POINTER :: cv(:,:), tv(:,:), tcv(:,:), tdiss(:,:)
  REAL(RFREAL), POINTER :: vol(:)

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ilev, iCOff, ijCOff, iNOff, ijNOff
  REAL(RFREAL), POINTER :: avgCo(:,:), sf(:,:), grad(:,:)
#endif
#ifdef RFLU
  INTEGER, POINTER      :: f2c(:,:)
  REAL(RFREAL), POINTER :: fn(:,:), grad(:,:,:)
  TYPE(t_patch), POINTER :: pPatch
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_RansSAVisFlux',&
  'TURB_RansSAVisFlux.F90' )

! get dimensions and pointers ------------------------------------------------

#ifdef RFLO
  ilev   =  region%currLevel
  cv     => region%levels(ilev)%mixt%cv
  tv     => region%levels(ilev)%mixt%tv
  tcv    => region%levels(ilev)%turb%cv
  tdiss  => region%levels(ilev)%turb%diss
  vol    => region%levels(iLev)%grid%vol
  
  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )
  iBegV = CV_SA_NUTIL
  iEndV = CV_SA_NUTIL
  iBegG = GR_SA_NUTILX
  iEndG = GR_SA_NUTILZ
#endif
#ifdef RFLU
  cv     => region%mixt%cv
  tv     => region%mixt%tv
  tcv    => region%turb%cv
  tdiss  => region%turb%diss  
  vol    => region%grid%vol
  iBegV = CV_SA_NUTIL
  iEndV = CV_SA_NUTIL
  iBegG = GR_SA_NUTILX
  iEndG = GR_SA_NUTILX
#endif

! get needed quantities

  cb2    = region%turbInput%const(MC_SA_CB2)
  rSigma = region%turbInput%const(MC_SA_RSIG)
  beta   = region%mixtInput%betrk(region%irkStep)*rSigma
  opcb2  = 1._RFREAL+cb2

! get gradients of tilde[nu]

#ifdef RFLO
  CALL RFLO_CalcGradVector( region,iBegV,iEndV,iBegG,iEndG, &
                                   region%levels(iLev)%turb%cv, &
                                   region%levels(iLev)%turb%gradi, &
                                   region%levels(iLev)%turb%gradj, &
                                   region%levels(iLev)%turb%gradk )
#endif
#ifdef RFLU
  CALL RFLU_ComputeGradFacesWrapper( region,iBegV,iEndV,iBegG,iEndG, &
                                     region%turb%cv,region%turb%gradi )
  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)
                                  
    IF ( RFLU_DecideNeedBGradFace(region,pPatch) .EQV. .TRUE. ) THEN
      CALL RFLU_ComputeGradBFacesWrapper( region,pPatch,iBegV,iEndV,iBegG,iEndG, &
                                          region%turb%cv,region%turb%bGradi )
    END IF ! RFLU_DecideNeedBGradFace
  END DO ! iPatch                                 
#endif

! interior fluxes -------------------------------------------------------------

  CALL ComputeFlux( DIRI )
#ifdef RFLO
  CALL ComputeFlux( DIRJ )
  CALL ComputeFlux( DIRK )
#endif

! fluxes through boundaries ---------------------------------------------------
#ifdef RFLO  
  DO ipatch=1,region%nPatches
    CALL TURB_RansSAVisFluxPatch( region,region%levels(ilev)%patches(ipatch) )
  ENDDO
#endif
#ifdef RFLU
  DO ipatch = 1,region%grid%nPatches
    CALL TURB_RansSAVisFluxPatch( region,region%patches(ipatch) )
  END DO ! ipatch
#endif

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

#ifdef RFLO
    IF (ijk==DIRI) THEN
      ibeg = ipcbeg+1
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = -1
      jadd = 0
      kadd = 0
      grad  => region%levels(ilev)%turb%gradi
      sf    => region%levels(ilev)%grid%si
      avgCo => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==DIRJ) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg+1
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = 0
      jadd = -1
      kadd = 0
      grad  => region%levels(ilev)%turb%gradj
      sf    => region%levels(ilev)%grid%sj
      avgCo => region%levels(iLev)%grid%c2fCoJ
    ELSEIF (ijk==DIRK) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg+1
      kend = kpcend
      iadd = 0
      jadd = 0
      kadd = -1
      grad  => region%levels(ilev)%turb%gradk
      sf    => region%levels(ilev)%grid%sk
      avgCo => region%levels(iLev)%grid%c2fCoK
    ENDIF
#endif
#ifdef RFLU
    ibeg =  1
    iend =  region%grid%nFaces
    f2c  => region%grid%f2c
    grad => region%turb%gradi
    fn   => region%grid%fn
#endif

! -- flux in ijk-direction (except through boundary) --------------------------

#ifdef RFLO
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
#endif
#ifdef RFLU
    ac0 = 0.5_RFREAL
    ac1 = 0.5_RFREAL
    DO ijkN = ibeg,iend
          ijkC0 = f2c(1,ijkN)
          ijkC1 = f2c(2,ijkN)
          sFace(1)= fn(XCOORD,ijkN)*fn(XYZMAG,ijkN)
          sFace(2)= fn(YCOORD,ijkN)*fn(XYZMAG,ijkN)
          sFace(3)= fn(ZCOORD,ijkN)*fn(XYZMAG,ijkN)           
#endif
          rhoa = ac0*cv(CV_MIXT_DENS ,ijkC0)+ac1*cv(CV_MIXT_DENS ,ijkC1)
          mua  = ac0*tv(TV_MIXT_MUEL ,ijkC0)+ac1*tv(TV_MIXT_MUEL ,ijkC1)
          rnua = ac0*tcv(CV_SA_NUTIL ,ijkC0)+ac1*tcv(CV_SA_NUTIL ,ijkC1)
          nuf  = (rnua + mua)/rhoa 
          nuc  = (tcv(CV_SA_NUTIL ,ijkC0) + tv(TV_MIXT_MUEL ,ijkC0))/ &
                  cv(CV_MIXT_DENS ,ijkC0) 

#ifdef RFLO
          nutilX = grad(GR_SA_NUTILX,ijkN)
          nutilY = grad(GR_SA_NUTILY,ijkN)
          nutilZ = grad(GR_SA_NUTILZ,ijkN)
#endif
#ifdef RFLU
          nutilX = grad(XCOORD,GR_SA_NUTILX,ijkN)
          nutilY = grad(YCOORD,GR_SA_NUTILX,ijkN)
          nutilZ = grad(ZCOORD,GR_SA_NUTILX,ijkN)
#endif

          fd = beta*(nuf*opcb2-nuc*cb2)* &
               (nutilX*sFace(1)+nutilY*sFace(2)+nutilZ*sFace(3))

          tdiss(CV_SA_NUTIL,ijkC0) = tdiss(CV_SA_NUTIL,ijkC0) + fd
          tdiss(CV_SA_NUTIL,ijkC1) = tdiss(CV_SA_NUTIL,ijkC1) - fd

#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
#else
    ENDDO       ! ijkN
#endif

  END SUBROUTINE ComputeFlux

END SUBROUTINE TURB_RansSAVisFlux

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_RansSAVisFlux.F90,v $
! Revision 1.12  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/08/19 15:40:38  mparmar
! Added use of RFLU_DecideNeedBGradFace
!
! Revision 1.9  2006/04/07 15:06:06  haselbac
! Bug fix: Incorrect ifs
!
! Revision 1.8  2006/04/07 14:56:02  haselbac
! Adapted to changes in f and bf grad routines
!
! Revision 1.7  2005/12/20 20:43:51  wasistho
! adapted to changing in Rocflu on face gradient routines
!
! Revision 1.6  2004/08/02 23:08:31  wasistho
! shift location of lines defining ac0 and ac1
!
! Revision 1.5  2004/08/02 21:55:39  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.4  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.3  2004/03/24 03:37:02  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.5  2003/10/25 22:07:26  wasistho
! modified non-conservative diffusion term
!
! Revision 1.4  2003/10/21 01:34:19  wasistho
! loop in computation of gradCell
!
! Revision 1.3  2003/10/20 20:28:21  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.2  2003/10/10 20:35:06  wasistho
! multiplied SA viscous fluxes by 1/sigma through beta
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







