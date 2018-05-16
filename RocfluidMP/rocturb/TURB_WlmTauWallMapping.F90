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
! Purpose: Mapping of modeled wall stress in body fitted coordinate to wall
!          stress components in Cartesian coordinate.
!
! Description: Stress involves second derivatives of velocities. The mapping
!              proceeds hence in two steps corresponding with inner and outer
!              derivatives. The body fitted coordinate system is patch%lbound
!              dependent. Eta axis always directs to region interior in wall
!              normal direction. Xi and zeta follows positive i, j, or k 
!              direction in the patch, and has the same directions for opposite
!              patches (lbound 1 and 2, 3 and 4, 5 and 6). Xi, eta, zeta is
!              thus right cycled at lbound 1, 3 and 5 and left cycled at 2, 4,
!              and 6. The face vector mapping is hence as follow:
!
!              lbound 1 : s_xi = -s_k, s_et = -si, s_zt = -s_j
!              lbound 2 : s_xi = -s_k, s_et = +si, s_zt = -s_j
!              lbound 3 : s_xi = -s_i, s_et = -sj, s_zt = -s_k
!              lbound 4 : s_xi = -s_i, s_et = +sj, s_zt = -s_k
!              lbound 5 : s_xi = -s_j, s_et = -sk, s_zt = -s_i
!              lbound 6 : s_xi = -s_j, s_et = +sk, s_zt = -s_i
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: modeled wall stresses in Cartesian coordinate.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_WlmTauWallMapping.F90,v 1.4 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmTauWallMapping( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices

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
  INTEGER :: i, j

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ilev, lbound, indxb, indxe, jndxb, jndxe
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag, ijkVal

  REAL(RFREAL), POINTER     :: vals(:,:)
  REAL(RFREAL), ALLOCATABLE :: tauWallTmp(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmTauWallMapping.F90,v $ $Revision: 1.4 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmTauWallMapping',&
  'TURB_WlmTauWallMapping.F90' )

! get dimensions and parameters -------------------------------------------------

  vals   => patch%valBola%vals 

#ifdef RFLO
  ilev   =  region%currLevel
  lbound =  patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )

! allocate temporary workspace and perform mapping ------------------------------

  n1    = ABS(patch%l1end-patch%l1beg)
  n2    = ABS(patch%l2end-patch%l2beg)
  iOff  = n1 + 1
  ijBeg = IndIJ( 0, 0,iOff)
  ijEnd = IndIJ(n1,n2,iOff)
#endif
#ifdef RFLU
  ijBeg =  1
  ijEnd =  patch%nBFaces
#endif

  ALLOCATE( tauWallTmp(ijBeg:ijEnd,TENSOR_ALL_NELM),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! get begin and end indices for i, j and k directions and check for consistency

#ifdef RFLO
  IF (lbound==1 .OR. lbound==2) THEN
    indxb = jbeg
    indxe = jend
    jndxb = kbeg
    jndxe = kend
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    indxb = kbeg
    indxe = kend
    jndxb = ibeg
    jndxe = iend
  ELSE
    indxb = ibeg
    indxe = iend
    jndxb = jbeg
    jndxe = jend
  ENDIF

  IF (n1/=(indxe-indxb) .OR. n2/=(jndxe-jndxb)) THEN
    CALL ErrorStop( global,ERR_PATCH_DIMENS,__LINE__, &
                   'Wlm patch dimension inconsistent' )
  ENDIF

  DO j=jndxb,jndxe
    DO i=indxb,indxe

      ijkVal = IndIJ(i-indxb ,j-jndxb ,indxe-indxb+1)
#endif
#ifdef RFLU
  DO i=ijBeg,ijEnd
      ijkVal = i
#endif

! --- First step:

      tauWallTmp(ijkVal,A11) = &
      vals(ijkVal,WLM_VALS_XIX)*vals(ijkVal,WLM_VALS_TAUUX) + &  
      vals(ijkVal,WLM_VALS_ETX)*vals(ijkVal,WLM_VALS_TAUUY) + &  
      vals(ijkVal,WLM_VALS_ZTX)*vals(ijkVal,WLM_VALS_TAUUZ)  

      tauWallTmp(ijkVal,A12) = &
      vals(ijkVal,WLM_VALS_XIY)*vals(ijkVal,WLM_VALS_TAUUX) + &  
      vals(ijkVal,WLM_VALS_ETY)*vals(ijkVal,WLM_VALS_TAUUY) + &  
      vals(ijkVal,WLM_VALS_ZTY)*vals(ijkVal,WLM_VALS_TAUUZ)  

      tauWallTmp(ijkVal,A13) = &
      vals(ijkVal,WLM_VALS_XIZ)*vals(ijkVal,WLM_VALS_TAUUX) + &  
      vals(ijkVal,WLM_VALS_ETZ)*vals(ijkVal,WLM_VALS_TAUUY) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*vals(ijkVal,WLM_VALS_TAUUZ)  
! -----
      tauWallTmp(ijkVal,A21) = &
      vals(ijkVal,WLM_VALS_XIX)*vals(ijkVal,WLM_VALS_TAUVX) + &  
      vals(ijkVal,WLM_VALS_ETX)*vals(ijkVal,WLM_VALS_TAUVY) + &  
      vals(ijkVal,WLM_VALS_ZTX)*vals(ijkVal,WLM_VALS_TAUVZ)  

      tauWallTmp(ijkVal,A22) = &
      vals(ijkVal,WLM_VALS_XIY)*vals(ijkVal,WLM_VALS_TAUVX) + &  
      vals(ijkVal,WLM_VALS_ETY)*vals(ijkVal,WLM_VALS_TAUVY) + &  
      vals(ijkVal,WLM_VALS_ZTY)*vals(ijkVal,WLM_VALS_TAUVZ)  

      tauWallTmp(ijkVal,A23) = &
      vals(ijkVal,WLM_VALS_XIZ)*vals(ijkVal,WLM_VALS_TAUVX) + &  
      vals(ijkVal,WLM_VALS_ETZ)*vals(ijkVal,WLM_VALS_TAUVY) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*vals(ijkVal,WLM_VALS_TAUVZ)  
! -----
      tauWallTmp(ijkVal,A31) = &
      vals(ijkVal,WLM_VALS_XIX)*vals(ijkVal,WLM_VALS_TAUWX) + &  
      vals(ijkVal,WLM_VALS_ETX)*vals(ijkVal,WLM_VALS_TAUWY) + &  
      vals(ijkVal,WLM_VALS_ZTX)*vals(ijkVal,WLM_VALS_TAUWZ)  

      tauWallTmp(ijkVal,A32) = &
      vals(ijkVal,WLM_VALS_XIY)*vals(ijkVal,WLM_VALS_TAUWX) + &  
      vals(ijkVal,WLM_VALS_ETY)*vals(ijkVal,WLM_VALS_TAUWY) + &  
      vals(ijkVal,WLM_VALS_ZTY)*vals(ijkVal,WLM_VALS_TAUWZ)  

      tauWallTmp(ijkVal,A33) = &
      vals(ijkVal,WLM_VALS_XIZ)*vals(ijkVal,WLM_VALS_TAUWX) + &  
      vals(ijkVal,WLM_VALS_ETZ)*vals(ijkVal,WLM_VALS_TAUWY) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*vals(ijkVal,WLM_VALS_TAUWZ)  

! --- Second step:

      vals(ijkVal,WLM_VALS_TAUUX) = &
      vals(ijkVal,WLM_VALS_XIX)*tauWallTmp(ijkVal,A11) + &  
      vals(ijkVal,WLM_VALS_ETX)*tauWallTmp(ijkVal,A21) + &  
      vals(ijkVal,WLM_VALS_ZTX)*tauWallTmp(ijkVal,A31)  

      vals(ijkVal,WLM_VALS_TAUUY) = &
      vals(ijkVal,WLM_VALS_XIX)*tauWallTmp(ijkVal,A12) + &  
      vals(ijkVal,WLM_VALS_ETX)*tauWallTmp(ijkVal,A22) + &  
      vals(ijkVal,WLM_VALS_ZTX)*tauWallTmp(ijkVal,A32)  

      vals(ijkVal,WLM_VALS_TAUUZ) = &
      vals(ijkVal,WLM_VALS_XIX)*tauWallTmp(ijkVal,A13) + &  
      vals(ijkVal,WLM_VALS_ETX)*tauWallTmp(ijkVal,A23) + &  
      vals(ijkVal,WLM_VALS_ZTX)*tauWallTmp(ijkVal,A33)  
! -----
      vals(ijkVal,WLM_VALS_TAUVX) = &
      vals(ijkVal,WLM_VALS_XIY)*tauWallTmp(ijkVal,A11) + &  
      vals(ijkVal,WLM_VALS_ETY)*tauWallTmp(ijkVal,A21) + &  
      vals(ijkVal,WLM_VALS_ZTY)*tauWallTmp(ijkVal,A31)  

      vals(ijkVal,WLM_VALS_TAUVY) = &
      vals(ijkVal,WLM_VALS_XIY)*tauWallTmp(ijkVal,A12) + &  
      vals(ijkVal,WLM_VALS_ETY)*tauWallTmp(ijkVal,A22) + &  
      vals(ijkVal,WLM_VALS_ZTY)*tauWallTmp(ijkVal,A32)  

      vals(ijkVal,WLM_VALS_TAUVZ) = &
      vals(ijkVal,WLM_VALS_XIY)*tauWallTmp(ijkVal,A13) + &  
      vals(ijkVal,WLM_VALS_ETY)*tauWallTmp(ijkVal,A23) + &  
      vals(ijkVal,WLM_VALS_ZTY)*tauWallTmp(ijkVal,A33)  
! -----
      vals(ijkVal,WLM_VALS_TAUWX) = &
      vals(ijkVal,WLM_VALS_XIZ)*tauWallTmp(ijkVal,A11) + &  
      vals(ijkVal,WLM_VALS_ETZ)*tauWallTmp(ijkVal,A21) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*tauWallTmp(ijkVal,A31)  

      vals(ijkVal,WLM_VALS_TAUWY) = &
      vals(ijkVal,WLM_VALS_XIZ)*tauWallTmp(ijkVal,A12) + &  
      vals(ijkVal,WLM_VALS_ETZ)*tauWallTmp(ijkVal,A22) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*tauWallTmp(ijkVal,A32)  

      vals(ijkVal,WLM_VALS_TAUWZ) = &
      vals(ijkVal,WLM_VALS_XIZ)*tauWallTmp(ijkVal,A13) + &  
      vals(ijkVal,WLM_VALS_ETZ)*tauWallTmp(ijkVal,A23) + &  
      vals(ijkVal,WLM_VALS_ZTZ)*tauWallTmp(ijkVal,A33)  

!wlmCheckprobe---------------------------------------------------------------
!      write(*,*) region%procId,patch%lbound,i,j, &
!                 vals(ijkVal,WLM_VALS_TAUUY),vals(ijkVal,WLM_VALS_TAUWY), &
!                 vals(ijkVal,WLM_VALS_TAUVX),vals(ijkVal,WLM_VALS_TAUVZ)
!----------------------------------------------------------------------------   
#ifdef RFLO
    ENDDO     ! i
  ENDDO       ! j
#endif
#ifdef RFLU
  ENDDO       ! i
#endif

! deallocate temporary workarrays

  DEALLOCATE( tauWallTmp )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmTauWallMapping

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmTauWallMapping.F90,v $
! Revision 1.4  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:01  wasistho
! changed nomenclature
!
! Revision 1.2  2004/03/02 03:51:04  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







