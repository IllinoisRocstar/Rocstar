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
! Purpose: Model wall stresses based on turbulent BL assumption with mixing
!          length model for eddy viscosity, and pressure gradient as well as
!          usteady terms taken into account.
!
! Description: Wall stresses are derived from modeling terms in BL equations.
!              Viscous term has been modeled in routine TURB_FloWlmUpdateLoglay
!              based on mixing lenth eddy viscosity model which equivalent to
!              log layer assumption at zero pressure gradient with surface
!              roughness included. BL convective term is neglected. Pressure 
!              gradient and time derivative term are computed here and added 
!              to the viscous model. The method is inspired by paper of
!              Hoffman and Benocci, "Approximate wall BC for LES", 5th
!              Advances in Turbulence, Siena, Italy 1994. Surface roughness
!              were not considered in their model.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: total and wall parallel components of wall stresses in body fitted 
!         coordinate.
!
! Notes: Unsteady term is not modeled yet. Its model term may be included 
!        in future.
!
!******************************************************************************
!
! $Id: TURB_WlmUpdateBndlay.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmUpdateBndlay( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#ifdef RFLO
#include "Indexing.h"
#endif

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: ijkVal

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ijBeg, ijEnd
  REAL(RFREAL)          :: wdist, dPdXi, dPdZt, utau, abVel
  REAL(RFREAL), POINTER :: vals(:,:)

#ifdef RFLO
  INTEGER :: n1, n2, iOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmUpdateBndlay.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmUpdateBndlay',&
  'TURB_WlmUpdateBndlay.F90' )

! get pointers ----------------------------------------------------------------

  vals  => patch%valBola%vals 

! get dimensions

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

! compute tau-wall and heat flux from total contributions

  DO ijkVal=ijBeg,ijEnd

    wdist = vals(ijkVal,WLM_VALS_WDIST)
    dPdXi = vals(ijkVal,WLM_VALS_DPDXI)
    dPdZt = vals(ijkVal,WLM_VALS_DPDZT)

    vals(ijkVal,WLM_VALS_TAUUY) = vals(ijkVal,WLM_VALS_TAUUY)-wdist*dPdXi     
    vals(ijkVal,WLM_VALS_TAUWY) = vals(ijkVal,WLM_VALS_TAUWY)-wdist*dPdZt

    vals(ijkVal,WLM_VALS_TAUVX) = vals(ijkVal,WLM_VALS_TAUUY)     
    vals(ijkVal,WLM_VALS_TAUVZ) = vals(ijkVal,WLM_VALS_TAUWY)      

! - store total wall stress in heat-flux array for heat transfer comp. later
    vals(ijkVal,WLM_VALS_HFLUX) = &
          SQRT( vals(ijkVal,WLM_VALS_TAUUY)*vals(ijkVal,WLM_VALS_TAUUY) + & 
                vals(ijkVal,WLM_VALS_TAUWY)*vals(ijkVal,WLM_VALS_TAUWY) )

! - Note: utau should not be derived from here (effect of pressure gradient 
!   exist), as it is used to initiate iteration on utau at the next stage 
!   in TURB_FloWlmUpdateLoglay (model for BL viscous term), in which pressure
!   gradient is not (yet) taken into account. On the contrary, for post-
!   processing/visualisation utau has to be extracted separately from wall 
!   stress tauWall computed in this routine: utau=sqrt(tauWall/rho)
   
!checkprobe--------------------------------------------------------------------
!    utau = SQRT( vals(ijkVal,WLM_VALS_HFLUX)/vals(ijkVal,WLM_VALS_DENS) )
!    abVel= SQRT( vals(ijkVal,WLM_VALS_XIV)**2+vals(ijkVal,WLM_VALS_ZTV)**2 )
!    write(*,*) region%procId,patch%lbound,ijkVal,vals(ijkVal,WLM_VALS_TAUUY), &
!               vals(ijkVal,WLM_VALS_TAUWY),utau/abVel,dPdXi,dPdZt
!------------------------------------------------------------------------------   

  ENDDO       ! ijkVal

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmUpdateBndlay

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmUpdateBndlay.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/24 03:37:03  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:01  wasistho
! changed nomenclature
!
! Revision 1.3  2004/03/02 03:51:27  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







