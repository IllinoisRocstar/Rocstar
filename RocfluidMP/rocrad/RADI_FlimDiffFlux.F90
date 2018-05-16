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
! Purpose: compute FLDTRAN diffusion flux: c*lamda(Er)/Kr*d_j(Er)
!
! Description: this routine compute d_j(Er), while the parameters defining
!              the diffusion coefficient are known from previous computations.
!
! Input: region  = data of current region
!
! Output: region%levels%radi%diss = diffusion flux added to FLD dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimDiffFlux.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimDiffFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, RFLO_CalcGradVector

#include "Indexing.h"
#endif
#ifdef RFLU
  USE RFLU_ModDifferentiation, ONLY: RFLU_ComputeGradFacesBak, & 
                                     RFLU_ComputeGradFacesPatchesBak
#endif
  USE RADI_ModInterfaces, ONLY : RADI_FlimDiffFluxPatch
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
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
  REAL(RFREAL) :: beta, sounda, flima, coefa, diffCoef, modSf, fd
  REAL(RFREAL) :: sFace(3), radEnX, radEnY, radEnZ
  REAL(RFREAL), POINTER :: dv(:,:), rcv(:,:), rdiss(:,:), flim(:), coef(:,:)
  REAL(RFREAL), POINTER :: qr(:)

#ifdef RFLO
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: ilev, iCOff, ijCOff, iNOff, ijNOff
  REAL(RFREAL), POINTER :: avgCo(:,:), sf(:,:), grad(:,:)
#endif
#ifdef RFLU
  INTEGER, POINTER      :: f2c(:,:)
  REAL(RFREAL), POINTER :: fn(:,:), grad(:,:,:)
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FlimDiffFlux',&
  'RADI_FlimDiffFlux.F90' )

! get dimensions and pointers ------------------------------------------------

#ifdef RFLO
  ilev  =  region%currLevel
  dv    => region%levels(ilev)%mixt%dv
  rcv   => region%levels(ilev)%radi%cv
  rdiss => region%levels(ilev)%radi%diss
  flim  => region%levels(iLev)%radi%fluxLim
  coef  => region%levels(iLev)%radi%radCoef
  
  CALL RFLO_GetDimensPhys( region,ilev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )
#endif
#ifdef RFLU
  dv    => region%mixt%dv
  rcv   => region%radi%cv
  rdiss => region%radi%diss  
  flim  => region%radi%fluxLim
  coef  => region%radi%radCoef
#endif
  iBegV = CV_RADI_ENER
  iEndV = CV_RADI_ENER
  iBegG = GR_RADI_EX
  iEndG = GR_RADI_EZ

! get needed quantities

  beta  = region%mixtInput%betrk(region%irkStep)

! get gradients of radiation Energy Er

#ifdef RFLO
  CALL RFLO_CalcGradVector( region,iBegV,iEndV,iBegG,iEndG, &
                                   region%levels(iLev)%radi%cv, &
                                   region%levels(iLev)%radi%gradi, &
                                   region%levels(iLev)%radi%gradj, &
                                   region%levels(iLev)%radi%gradk )
#endif
#ifdef RFLU
  CALL RFLU_ComputeGradFacesBak( region,iBegV,iEndV,iBegG,iEndG, &
                                 region%radi%cv,region%radi%gradi )
  CALL RFLU_ComputeGradFacesPatchesBak( region,iBegV,iEndV,iBegG,iEndG, &
                                        region%radi%cv,region%radi%bGradi )
#endif

! interior fluxes -------------------------------------------------------------

  CALL ComputeFlux( ICOORD )
#ifdef RFLO
  CALL ComputeFlux( JCOORD )
  CALL ComputeFlux( KCOORD )
#endif

! fluxes through boundaries ---------------------------------------------------
#ifdef RFLO  
  DO ipatch=1,region%nPatches
    CALL RADI_FlimDiffFluxPatch( region,region%levels(ilev)%patches(ipatch) )
  ENDDO
#endif
#ifdef RFLU
  DO ipatch = 1,region%grid%nPatches
    CALL RADI_FlimDiffFluxPatch( region,region%patches(ipatch) )
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
    IF (ijk==ICOORD) THEN
      ibeg = ipcbeg+1
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = -1
      jadd = 0
      kadd = 0
      qr    => region%levels(iLev)%radi%qri
      grad  => region%levels(ilev)%radi%gradi
      sf    => region%levels(ilev)%grid%si
      avgCo => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==JCOORD) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg+1
      jend = jpcend
      kbeg = kpcbeg
      kend = kpcend
      iadd = 0
      jadd = -1
      kadd = 0
      qr    => region%levels(iLev)%radi%qrj
      grad  => region%levels(ilev)%radi%gradj
      sf    => region%levels(ilev)%grid%sj
      avgCo => region%levels(iLev)%grid%c2fCoJ
    ELSEIF (ijk==KCOORD) THEN
      ibeg = ipcbeg
      iend = ipcend
      jbeg = jpcbeg
      jend = jpcend
      kbeg = kpcbeg+1
      kend = kpcend
      iadd = 0
      jadd = 0
      kadd = -1
      qr    => region%levels(iLev)%radi%qrk
      grad  => region%levels(ilev)%radi%gradk
      sf    => region%levels(ilev)%grid%sk
      avgCo => region%levels(iLev)%grid%c2fCoK
    ENDIF
#endif
#ifdef RFLU
    ibeg =  1
    iend =  region%grid%nFaces
    f2c  => region%grid%f2c
    qr   => region%radi%qri
    grad => region%radi%gradi
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
          sounda   = ac0*dv(DV_MIXT_SOUN ,ijkC0)+ac1*dv(DV_MIXT_SOUN ,ijkC1)
          flima    = ac0*flim(ijkC0)+ac1*flim(ijkC1)
          coefa    = ac0*coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                     ac1*coef(ijkC1,RADI_COEFF_EXTINCT)
          diffCoef = sounda*flima/coefa 

#ifdef RFLO
          radEnX = grad(GR_RADI_EX,ijkN)
          radEnY = grad(GR_RADI_EY,ijkN)
          radEnZ = grad(GR_RADI_EZ,ijkN)
#endif
#ifdef RFLU
          radEnX = grad(XCOORD,GR_RADI_EX,ijkN)
          radEnY = grad(YCOORD,GR_RADI_EX,ijkN)
          radEnZ = grad(ZCOORD,GR_RADI_EX,ijkN)
#endif

          fd = diffCoef* &
               (radEnX*sFace(1)+radEnY*sFace(2)+radEnZ*sFace(3))

          rdiss(CV_RADI_ENER,ijkC0) = rdiss(CV_RADI_ENER,ijkC0) + fd*beta
          rdiss(CV_RADI_ENER,ijkC1) = rdiss(CV_RADI_ENER,ijkC1) - fd*beta

! ------- store radiant flux normal to cell face in positive direction in qr
          modSf     = SQRT( sFace(1)*sFace(1) + &
                            sFace(2)*sFace(2) + &
                            sFace(3)*sFace(3))
          qr(ijkN) = fd/modSf
#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
#else
    ENDDO       ! ijkN
#endif

  END SUBROUTINE ComputeFlux

END SUBROUTINE RADI_FlimDiffFlux

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RADI_FlimDiffFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







