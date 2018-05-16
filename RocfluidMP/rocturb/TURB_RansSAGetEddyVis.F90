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
! Purpose: update eddy viscosity and turbulent thermal conductivity from 
!          working variable tilde[nu]
!
! Description: in SA formulation (La Recherche Aerospatiale 1, 5 (1994)):
!              xi   = tilde[nu]/nu
!              fv1  = xi^3/(xi^3 + cv1^3)
!              nu_t = tilde[nu]*fv1
!
! Input: region = data of current region,
!
! Output: region%levels%mixt%tv(MUET,:) and tv(TCOT,:)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansSAGetEddyVis.F90,v 1.6 2009/08/12 04:15:59 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansSAGetEddyVis( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic
  
! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: indCp, ibc, iec, npow
  REAL(RFREAL), POINTER :: cv(:,:), tv(:,:), gv(:,:), tcv(:,:)
  REAL(RFREAL)          :: rPrt, cpPrt, xi, fv1, cv1, rnuet

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_RansSAGetEddyVis',&
  'TURB_RansSAGetEddyVis.F90' )

! get coefficients and constants ----------------------------------------------

  cv1 = region%turbInput%const(MC_SA_CV1)
  
  IF (region%turbInput%functV1 == SA_FV1_POW3) THEN
    npow = 3
  ELSE
    npow = 2
  ENDIF

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel
  indCp = region%levels(iLev)%mixt%indCp
  rPrt  = 1._RFREAL/region%levels(iLev)%mixt%prTurb

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv      => region%levels(iLev)%mixt%cv
  tv      => region%levels(iLev)%mixt%tv
  gv      => region%levels(iLev)%mixt%gv
  tcv     => region%levels(iLev)%turb%cv
#endif

#ifdef RFLU
  indCp = region%mixtInput%indCp
  rPrt  = 1._RFREAL/region%mixtInput%prTurb

  ibc = 1
  iec = region%grid%nCellsTot        

  cv      => region%mixt%cv
  tv      => region%mixt%tv
  gv      => region%mixt%gv
  tcv     => region%turb%cv
#endif

! update eddy viscosity and turbulent thermal conductivity

  DO ic=ibc,iec
    cpPrt = gv(GV_MIXT_CP,ic*indCp)*rPrt

! Temporary clipping fix
    tcv(CV_SA_NUTIL,ic) = MAX(tcv(CV_SA_NUTIL,ic),REAL_SMALL)
!
    rnuet = tcv(CV_SA_NUTIL,ic)
    xi    = rnuet/tv(TV_MIXT_MUEL,ic)
    fv1   = xi**npow/(xi**npow + cv1**npow)
    tv(TV_MIXT_MUET,ic) = fv1*rnuet
    tv(TV_MIXT_TCOT,ic) = cpPrt*tv(TV_MIXT_MUET,ic)
  ENDDO !ic

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansSAGetEddyVis

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansSAGetEddyVis.F90,v $
! Revision 1.6  2009/08/12 04:15:59  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.5  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
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
! Revision 1.3  2004/02/19 04:04:05  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.2  2003/10/20 20:27:56  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







