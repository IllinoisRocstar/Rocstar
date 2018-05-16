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
!          through a patch
!
! Description: this routine works in the same way as RADI_FlimDiffFlux
!              but applied on region patches
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: region%levels%radi%diss = diffusion flux added to FLD dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimDiffFluxPatch.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimDiffFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, iC

! ... local variables
  INTEGER :: bcType, ijkCB0, ijkCD, ijkNB
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL) :: beta, sounda, flima, coefa, diffCoef, modSf
  REAL(RFREAL) :: radEnX, radEnY, radEnZ, fd, sf(3), ac0, ac1
  REAL(RFREAL), POINTER :: dv(:,:), rcv(:,:), rdiss(:,:), flim(:), coef(:,:)
  REAL(RFREAL), POINTER :: qr(:)

#ifdef RFLO
  INTEGER :: inode, jnode, knode, idir, jdir, kdir
  INTEGER :: ilev, lbound, iCOff, ijCOff, iNOff, ijNOff, acId0, acId1
  REAL(RFREAL)          :: sgn
  REAL(RFREAL), POINTER :: avgCo(:,:), sFace(:,:), grad(:,:)
#endif
#ifdef RFLU
  INTEGER :: ifgBeg, ijkNBG
  REAL(RFREAL), POINTER :: fn(:,:), grad(:,:,:) 
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FlimDiffFluxPatch',&
  'RADI_FlimDiffFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  bcType = patch%bcType

#ifdef RFLO
  ilev   = region%currLevel
  lbound = patch%lbound
  
  dv    => region%levels(ilev)%mixt%dv
  rcv   => region%levels(ilev)%radi%cv
  rdiss => region%levels(ilev)%radi%diss  
  flim  => region%levels(iLev)%radi%fluxLim
  coef  => region%levels(iLev)%radi%radCoef
#endif  
#ifdef RFLU
  dv    => region%mixt%dv
  rcv   => region%radi%cv
  rdiss => region%radi%diss   
  flim  => region%radi%fluxLim
  coef  => region%radi%radCoef
#endif

! get coefficients -----------------------------------------------------------

  beta   = region%mixtInput%betrk(region%irkStep)

#ifdef RFLO

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
    grad  => region%levels(ilev)%radi%gradi
    qr    => region%levels(ilev)%radi%qri
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(ilev)%grid%sj
    grad  => region%levels(ilev)%radi%gradj
    qr    => region%levels(ilev)%radi%qrj
  ELSE
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(ilev)%grid%sk
    grad  => region%levels(ilev)%radi%gradk
    qr    => region%levels(ilev)%radi%qrk
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
          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
#endif
#ifdef RFLU
    ibeg   = 1
    iend   = patch%nBFaces
! TEMPORARY : removing usage of bf2bg from everywhere
!    ifgBeg = patch%bf2bg(BF2BG_BEG)
    ac0    = 0.5_RFREAL
    ac1    = 0.5_RFREAL

    grad => region%radi%bGradi
    fn   => patch%fn

    DO iC=ibeg,iend
          ijkCB0 = patch%bf2c(iC)
          ijkCD  = ijkCB0
          ijkNB  = iC
          ijkNBG = iC + ifgBeg-1

          sf(1)  = fn(XCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(2)  = fn(YCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(3)  = fn(ZCOORD,ijkNB)*fn(XYZMAG,ijkNB)             
#endif
          sounda   = ac0*dv(DV_MIXT_SOUN,ijkCB0)+ac1*dv(DV_MIXT_SOUN,ijkCD)
          flima    = ac0*flim(ijkCB0)+ac1*flim(ijkCD)
          coefa    = ac0*coef(ijkCB0,RADI_COEFF_EXTINCT)+ &
                     ac1*coef(ijkCD,RADI_COEFF_EXTINCT)
          diffCoef = sounda*flima/coefa 

#ifdef RFLO
          radEnX = grad(GR_RADI_EX,ijkNB)
          radEnY = grad(GR_RADI_EY,ijkNB)
          radEnZ = grad(GR_RADI_EZ,ijkNB)
#endif
#ifdef RFLU
          radEnX = grad(XCOORD,GR_RADI_EX,ijkNBG)
          radEnY = grad(YCOORD,GR_RADI_EX,ijkNBG)
          radEnZ = grad(ZCOORD,GR_RADI_EX,ijkNBG)
#endif

          fd = diffCoef* &
               (radEnX*sf(1)+radEnY*sf(2)+radEnZ*sf(3))

          rdiss(CV_RADI_ENER,ijkCB0) = rdiss(CV_RADI_ENER,ijkCB0)+fd*beta

! ------- store rad. flux normal to surf. in positive direction (hence sgn* )
          modSf     = SQRT( sf(1)*sf(1) + sf(2)*sf(2) + sf(3)*sf(3) )  
          qr(ijkNB) = sgn*fd/modSf
#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
  ENDIF
#endif 
      
#ifdef RFLU
  ENDDO         ! iC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FlimDiffFluxPatch

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RADI_FlimDiffFluxPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:15  mparmar
! Removed bf2bg
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







