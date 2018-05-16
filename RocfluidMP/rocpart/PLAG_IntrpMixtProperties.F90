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
! Purpose: Interpolate mixture properties onto particle locations
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_IntrpMixtProperties.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_intrpMixtProperties( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset 
#endif
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#ifdef RFLO
#include "Indexing.h"
#endif

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  
! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: ibc,iCOff,idcbeg,idcend,iec,ijCOff,iLev,   &
             ipcbeg,ipcend,jdcbeg,jdcend,jpcbeg,jpcend, &
             kdcbeg,kdcend,kpcbeg,kpcend
#endif 
  INTEGER :: iCell, nPcls
  INTEGER, POINTER, DIMENSION(:,:)    :: pAiv

  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv, pDv, pTv
#ifdef RFLO  
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCofg, pCvMixt, pDvMixt, pTvMixt
#endif
#ifdef RFLU  
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCvMixt, pDvMixt, pTvMixt
#endif
  
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 
 
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_IntrpMixtProperties.F90,v $ $Revision: 1.4 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_IntrpMixtProperties',&
  'PLAG_IntrpMixtProperties.F90' )

! Get dimensions --------------------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel
 
  nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  nPcls = region%plag%nPcls
#endif

! Set pointers ----------------------------------------------------------------

#ifdef RFLO
  pMixt => region%levels(iLev)%mixt
  pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pMixt => region%mixt
  pPlag => region%plag
#endif

! Set pointers for mixture ----------------------------------------------------

#ifdef RFLO  
  pCofg   => region%levels(iLev)%grid%cofg 
#endif

  pCvMixt => pMixt%cv
  pDvMixt => pMixt%dv
  pTvMixt => pMixt%tv

! Set pointers for discrete particles -----------------------------------------

  pAiv    => pPlag%aiv
  pCv     => pPlag%cv  
  pDv     => pPlag%dv
  pTv     => pPlag%tv

#ifdef RFLO     
! Get cell dimensions ---------------------------------------------------------
  
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset(  region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
       
! Get physical dimensions  ----------------------------------------------------
      
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
#endif

#ifdef RFLU
! Check that have primitive state vector --------------------------------------

  IF ( pMixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region%mixt%cvState
#endif  
     
! Assume a piece-wise constant distribution in a cell -------------------------
!  Note: Distribution has been properly corrected for RFLU
!        with the routine PLAG_RFLU_CorrectMixtProperties.
!        This routine is called after the cell gradients are computed
    
  DO iPcls = 1, nPcls

    iCell = pAiv(AIV_PLAG_ICELLS,iPcls)
#ifdef RFLO
    pDv(DV_PLAG_UVELMIXT,iPcls) = pDvMixt(DV_MIXT_UVEL,iCell)
    pDv(DV_PLAG_VVELMIXT,iPcls) = pDvMixt(DV_MIXT_VVEL,iCell)  
    pDv(DV_PLAG_WVELMIXT,iPcls) = pDvMixt(DV_MIXT_WVEL,iCell)
    pDv(DV_PLAG_PRESMIXT,iPcls) = pDvMixt(DV_MIXT_PRES,iCell)
    pDv(DV_PLAG_TEMPMIXT,iPcls) = pDvMixt(DV_MIXT_TEMP,iCell)
    pDv(DV_PLAG_DENSMIXT,iPcls) = pCvMixt(CV_MIXT_DENS,iCell)
#endif
#ifdef RFLU
    pDv(DV_PLAG_UVELMIXT,iPcls) = pCvMixt(CV_MIXT_XVEL,iCell)
    pDv(DV_PLAG_VVELMIXT,iPcls) = pCvMixt(CV_MIXT_YVEL,iCell)  
    pDv(DV_PLAG_WVELMIXT,iPcls) = pCvMixt(CV_MIXT_ZVEL,iCell)
    pDv(DV_PLAG_PRESMIXT,iPcls) = pDvMixt(DV_MIXT_PRES,iCell)
    pDv(DV_PLAG_TEMPMIXT,iPcls) = pDvMixt(DV_MIXT_TEMP,iCell)
    pDv(DV_PLAG_DENSMIXT,iPcls) = pCvMixt(CV_MIXT_DENS,iCell)
#endif

    pTv(TV_PLAG_MUELMIXT,iPcls) = pTvMixt(TV_MIXT_MUEL,iCell)
    pTv(TV_PLAG_TCOLMIXT,iPcls) = pTvMixt(TV_MIXT_TCOL,iCell)
  END DO ! iPcls
   
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_intrpMixtProperties

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_IntrpMixtProperties.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/11/30 22:18:57  fnajjar
! Changed bcos of corr routine, now only do piecewise const here
!
! Revision 1.1  2004/12/01 20:57:51  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/02/26 21:02:15  haselbac
! Fixed RFLU bug, removed superfluous statements
!
! Revision 1.4  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/01/16 20:18:53  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







