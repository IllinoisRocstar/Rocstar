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
!          through a patch
!
! Description: this routine works in the same way as TURB_RansSAVisFlux
!              but applied on region patches
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: region%levels%turb%diss = viscous flux added to SA dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansSAVisFluxPatch.F90,v 1.10 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansSAVisFluxPatch( region,patch )

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
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, iC

! ... local variables
  INTEGER :: bcType, ijkCB0, ijkCD, ijkNB
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL)          :: cb2, opcb2, rSigma, beta, rhoa, mua, rnua, nuf, nuc
  REAL(RFREAL)          :: nutilX, nutilY, nutilZ, fd, sf(3), ac0, ac1
  REAL(RFREAL), POINTER :: cv(:,:), tv(:,:), tcv(:,:), tdiss(:,:)

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

  CALL RegisterFunction( region%global,'TURB_RansSAVisFluxPatch',&
  'TURB_RansSAVisFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  bcType = patch%bcType

#ifdef RFLO
  ilev   = region%currLevel
  lbound = patch%lbound
  
  cv   => region%levels(ilev)%mixt%cv
  tv   => region%levels(ilev)%mixt%tv
  tcv  => region%levels(ilev)%turb%cv
  tdiss=> region%levels(ilev)%turb%diss  
#endif  
#ifdef RFLU
  cv   => region%mixt%cv
  tv   => region%mixt%tv
  tcv  => region%turb%cv
  tdiss=> region%turb%diss   
#endif

! get coefficients -----------------------------------------------------------

  cb2    = region%turbInput%const(MC_SA_CB2)
  rSigma = region%turbInput%const(MC_SA_RSIG)
  beta   = region%mixtInput%betrk(region%irkStep)*rSigma
  opcb2  = 1._RFREAL+cb2

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
    grad  => region%levels(ilev)%turb%gradi
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    sFace => region%levels(ilev)%grid%sj
    grad  => region%levels(ilev)%turb%gradj
  ELSE
    avgCo => region%levels(iLev)%grid%c2fCoK
    sFace => region%levels(ilev)%grid%sk
    grad  => region%levels(ilev)%turb%gradk
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

    grad => region%turb%bGradi
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
          rhoa = ac0*cv(CV_MIXT_DENS ,ijkCB0)+ac1*cv(CV_MIXT_DENS ,ijkCD)
          mua  = ac0*tv(TV_MIXT_MUEL ,ijkCB0)+ac1*tv(TV_MIXT_MUEL ,ijkCD)
          rnua = ac0*tcv(CV_SA_NUTIL ,ijkCB0)+ac1*tcv(CV_SA_NUTIL ,ijkCD)
          nuf  = (rnua + mua)/rhoa
          nuc  = (tcv(CV_SA_NUTIL ,ijkCB0) + tv(TV_MIXT_MUEL ,ijkCB0))/ &
                  cv(CV_MIXT_DENS ,ijkCB0) 

#ifdef RFLO
          nutilX = grad(GR_SA_NUTILX,ijkNB)
          nutilY = grad(GR_SA_NUTILY,ijkNB)
          nutilZ = grad(GR_SA_NUTILZ,ijkNB)
#endif
#ifdef RFLU
          nutilX = grad(XCOORD,GR_SA_NUTILX,ijkNBG)
          nutilY = grad(YCOORD,GR_SA_NUTILX,ijkNBG)
          nutilZ = grad(ZCOORD,GR_SA_NUTILX,ijkNBG)
#endif

          fd = beta*(nuf*opcb2-nuc*cb2)* &
               (nutilX*sf(1)+nutilY*sf(2)+nutilZ*sf(3))

          tdiss(CV_SA_NUTIL,ijkCB0) = tdiss(CV_SA_NUTIL,ijkCB0)+fd

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

END SUBROUTINE TURB_RansSAVisFluxPatch

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_RansSAVisFluxPatch.F90,v $
! Revision 1.10  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/08/19 15:41:01  mparmar
! Removed bf2bg
!
! Revision 1.7  2004/08/04 22:07:46  wasistho
! bugfixed: ac0 and ac1 are common to flo and flu
!
! Revision 1.6  2004/08/02 23:08:35  wasistho
! shift location of lines defining ac0 and ac1
!
! Revision 1.5  2004/08/02 21:55:46  wasistho
! replaced cell2face midpoint by linear averaging
!
! Revision 1.4  2004/03/25 04:40:41  wasistho
! prepared for RFLU
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
! Revision 1.4  2003/10/25 22:07:34  wasistho
! modified non-conservative diffusion term
!
! Revision 1.3  2003/10/20 20:28:28  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.2  2003/10/10 20:35:13  wasistho
! multiplied SA viscous fluxes by 1/sigma through beta
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







