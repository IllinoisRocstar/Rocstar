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
! Purpose: Set initial values of variables needed in wlm computations.
!
! Description: Friction velocity u_tau is set to 0.05 u_edge.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: initial u_tau.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_WlmInitia.F90,v 1.5 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmInitia( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset

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

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ijkC, ijkVal
  REAL(RFREAL)          :: u, v, w
  REAL(RFREAL), POINTER :: vals(:,:)

#ifdef RFLO
  INTEGER :: ilev, lbound, iCOff, ijCOff
  REAL(RFREAL), POINTER :: dv(:,:)
#endif
#ifdef RFLU
  REAL(RFREAL)          :: rDens
  REAL(RFREAL), POINTER :: cv(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmInitia.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmInitia',&
  'TURB_WlmInitia.F90' )

! get common pointers ---------------------------------------------------------

  vals => patch%valBola%vals 

#ifdef RFLO
! get dimensions and parameters -----------------------------------------------

  ilev   =  region%currLevel
  lbound =  patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

! get pointers

  dv   => region%levels(ilev)%mixt%dv

! estimate initial friction velocity utau

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend

        ijkC = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  

        IF (lbound==1 .OR. lbound==2) THEN
          ijkVal = IndIJ(j-jbeg ,k-kbeg ,jend-jbeg+1)
        ELSEIF (lbound==3 .OR. lbound==4) THEN
          ijkVal = IndIJ(k-kbeg ,i-ibeg ,kend-kbeg+1)
        ELSEIF (lbound==5 .OR. lbound==6) THEN
          ijkVal = IndIJ(i-ibeg ,j-jbeg ,iend-ibeg+1)
        ENDIF

        u = dv(DV_MIXT_UVEL,ijkC)     
        v = dv(DV_MIXT_VVEL,ijkC)     
        w = dv(DV_MIXT_WVEL,ijkC)
#endif
#ifdef RFLU
! specific Rocflu, check the state of cv first --------------------------------
  IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) &
                            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)

! get dimensions and parameters -----------------------------------------------

  ibeg =  1
  iend =  patch%nBFaces

! get pointers
  
  cv   => region%mixt%cv 

! estimate initial friction velocity utau

  DO i=ibeg,iend
        ijkC   = patch%bf2c(i) 
        ijkVal = i

        rDens  = 1._RFREAL/cv(CV_MIXT_DENS,ijkC) 
        u = cv(CV_MIXT_XMOM,ijkC)*rDens
        v = cv(CV_MIXT_YMOM,ijkC)*rDens
        w = cv(CV_MIXT_ZMOM,ijkC)*rDens
#endif
        vals(ijkVal,WLM_VALS_UTAU) = 0.05*SQRT(u*u+v*v+w*w)

#ifdef RFLO
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  ENDDO       ! i
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmInitia

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmInitia.F90,v $
! Revision 1.5  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/24 03:37:03  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:01  wasistho
! changed nomenclature
!
! Revision 1.2  2004/03/02 03:49:45  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







