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
! Purpose: Compute total viscous and heat fluxes through current patch using 
!          wall stress and heat transfer from wall layer model selected.
!
! Description: This routine works in the same way as TURB_VisFluxEddyPatch
!              or TURB_VFluxHybridPatch but using wall stress from wall layer 
!              model instead of mu*Sij or SGS-tau_ij, respectively.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: mixt%diss = total viscous fluxes added to dissipation
!
! Notes: Unlike the general viscous flux patch routines (TURB_visFluxEdduPatch
!        and TURB_vFluxHybridPatch), this wlm flux routine coverts both eddy
!        viscosity model (EVM, fixed and dynamic Smagorinsky) as well as
!        turbulent stress models (TSM, scale similarity and dynamic mixed).
!
!******************************************************************************
!
! $Id: TURB_WlmFluxPatch.F90,v 1.4 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices,RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset,RFLO_GetNodeOffset
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, iC

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: bcType, ijkCB0, ijkNB, ijkVal
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  REAL(RFREAL)          :: beta, heatFlux
  REAL(RFREAL)          :: fd(4),sf(3),twij(3,3)
  REAL(RFREAL), POINTER :: diss(:,:),vals(:,:)

#ifdef RFLO
  INTEGER :: inode, jnode, knode, idir, jdir, kdir
  INTEGER :: ilev, lbound, iCOff, ijCOff, iNOff, ijNOff
  REAL(RFREAL)          :: sgn
  REAL(RFREAL), POINTER :: sFace(:,:)
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER :: fn(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmFluxPatch.F90,v $ $Revision: 1.4 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmFluxPatch',&
  'TURB_WlmFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  vals   => patch%valBola%vals

#ifdef RFLO
  ilev   =  region%currLevel
  lbound =  patch%lbound
  diss   => region%levels(ilev)%mixt%diss  
#endif  
#ifdef RFLU
  diss   => region%mixt%diss  
#endif  

! get coefficients -----------------------------------------------------------

  bcType = patch%bcType
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
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF

! get the appropriate face vector --------------------------------------------

  IF (lbound==1 .OR. lbound==2) THEN
    sFace => region%levels(ilev)%grid%si
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sFace => region%levels(ilev)%grid%sj
  ELSE
    sFace => region%levels(ilev)%grid%sk
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
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)  ! bnd nodes
          ijkVal = ABS(idir)*IndIJ(j-jbeg,k-kbeg,jend-jbeg+1) + &
                   ABS(jdir)*IndIJ(k-kbeg,i-ibeg,kend-kbeg+1) + &
                   ABS(kdir)*IndIJ(i-ibeg,j-jbeg,iend-ibeg+1)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)
#endif
#ifdef RFLU
    ibeg  =  1
    iend  =  patch%nBFaces
    fn    => patch%fn

    DO iC=ibeg,iend
          ijkCB0 = patch%bf2c(iC)
          ijkNB  = iC
          ijkVal = iC

          sf(1)  = fn(XCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(2)  = fn(YCOORD,ijkNB)*fn(XYZMAG,ijkNB)
          sf(3)  = fn(ZCOORD,ijkNB)*fn(XYZMAG,ijkNB)           
#endif

          twij(1,1) = vals(ijkVal,WLM_VALS_TAUUX)
          twij(1,2) = vals(ijkVal,WLM_VALS_TAUUY)
          twij(1,3) = vals(ijkVal,WLM_VALS_TAUUZ)

          twij(2,1) = vals(ijkVal,WLM_VALS_TAUVX)
          twij(2,2) = vals(ijkVal,WLM_VALS_TAUVY)
          twij(2,3) = vals(ijkVal,WLM_VALS_TAUVZ)

          twij(3,1) = vals(ijkVal,WLM_VALS_TAUWX)
          twij(3,2) = vals(ijkVal,WLM_VALS_TAUWY)
          twij(3,3) = vals(ijkVal,WLM_VALS_TAUWZ)

          heatFlux  = SQRT(sf(1)*sf(1)+sf(2)*sf(2)+sf(3)*sf(3))* &
                      vals(ijkVal,WLM_VALS_HFLUX) 

          fd(1) = twij(1,1)*sf(1)+twij(1,2)*sf(2)+twij(1,3)*sf(3)
          fd(2) = twij(2,1)*sf(1)+twij(2,2)*sf(2)+twij(2,3)*sf(3)
          fd(3) = twij(3,1)*sf(1)+twij(3,2)*sf(2)+twij(3,3)*sf(3)
          fd(4) = heatFlux

!wlmCheckprobe---------------------------------------------------------------
!    write(*,*)region%procId,i,j,k,twij(1,2),twij(3,2),twij(2,1),twij(2,3), &
!              heatFlux
!----------------------------------------------------------------------------

          diss(CV_MIXT_XMOM,ijkCB0) = diss(CV_MIXT_XMOM,ijkCB0)+fd(1)*beta
          diss(CV_MIXT_YMOM,ijkCB0) = diss(CV_MIXT_YMOM,ijkCB0)+fd(2)*beta
          diss(CV_MIXT_ZMOM,ijkCB0) = diss(CV_MIXT_ZMOM,ijkCB0)+fd(3)*beta
          diss(CV_MIXT_ENER,ijkCB0) = diss(CV_MIXT_ENER,ijkCB0)+fd(4)*beta

#ifdef RFLO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
  ENDIF         ! bcType
#endif
#ifdef RFLU
  ENDDO         ! iC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmFluxPatch

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmFluxPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/24 03:37:03  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2004/03/02 03:49:30  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







