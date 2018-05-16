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
! Purpose: add laminar and turbulent transport variables to form total tv
!
! Description: in this routine we compute
!              mu_tot  = mu_l  + mu_t  with l=laminar, t=turbulent
!              tco_tot = tco_l + tco_t  with l=laminar, t=turbulent
!
! Input: region  = data of current region
!        indxMu  = component index for total dynamic viscosity Mu 
!        indxTCo = component index for total thermal conductivity TCo 
!        tvt     = transport variables, mu and tco
!
! Output: region%levels%turb%tv = total tv
!
! Notes: none
!
!******************************************************************************
!
! $Id: TURB_RansTotalTv.F90,v 1.4 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansTotalTv( region,indxMu,indxTCo,tvt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: indxMu, indxTCo
  REAL(RFREAL), POINTER :: tvt(:,:)

! ... loop variables
  INTEGER :: iC

! ... local variables
  INTEGER :: ibc, iec
  REAL(RFREAL), POINTER :: tv(:,:)

#ifdef RFLO
  INTEGER :: iLev, iCOff,ijCOff, idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_RansTotalTv',&
  'TURB_RansTotalTv.F90' )

! get dimensions and pointers ------------------------------------------------

#ifdef RFLO
  iLev =  region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  tv  => region%levels(ilev)%mixt%tv
#endif
#ifdef RFLU
  ibc =  1
  iec =  region%grid%nCellsTot
  tv  => region%mixt%tv
#endif

  DO iC=ibc,iec
    tvt(indxMu ,iC) = tv(TV_MIXT_MUEL,iC)+tv(TV_MIXT_MUET,iC)
    tvt(indxTCo,iC) = tv(TV_MIXT_TCOL,iC)+tv(TV_MIXT_TCOT,iC)
  ENDDO ! iC

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_RansTotalTv

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_RansTotalTv.F90,v $
! Revision 1.4  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2003/10/08 03:55:27  wasistho
! bug fixed MUET to TCOL
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







