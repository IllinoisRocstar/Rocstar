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
! Purpose: added source terms of SA/DES to the rhs.
!
! Description: see formulation of SA (La Recherche Aerospatiale 1, 5 (1994))
!
! Input: region = data of current region,
!
! Output: region%levels%turb%rhs , source terms added to rhs
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansSASourceTerms.F90,v 1.5 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansSASourceTerms( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  USE TURB_ModInterfaces, ONLY : TURB_CalcVortic
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ijkC
  
! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: ibc, iec
  REAL(RFREAL), POINTER :: cv(:,:), tv(:,:), vol(:), wdist(:)
  REAL(RFREAL), POINTER :: tcv(:,:), trhs(:,:), vort(:,:), dsterm(:,:)
  REAL(RFREAL) :: one6th, cv1, cw1, cw2, cw3, cb1, rKappa, rSigma
  REAL(RFREAL) :: rnuet, nuet, xi, fv1, fv2, rWdist, sTilde, ro, go, fw
  REAL(RFREAL) :: vortMag, prod, destr, sTerm

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff,ijCOff
#endif

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'TURB_RansSASourceTerms',&
  'TURB_RansSASourceTerms.F90' )

! get coefficients and parameters ---------------------------------------------

  one6th = 1._RFREAL/6._RFREAL  
  cv1    = region%turbInput%const(MC_SA_CV1)
  cw1    = region%turbInput%const(MC_SA_CW1)
  cw2    = region%turbInput%const(MC_SA_CW2)
  cw3    = region%turbInput%const(MC_SA_CW3)
  cb1    = region%turbInput%const(MC_SA_CB1)
  rKappa = region%turbInput%const(MC_SA_RKAP)
  rSigma = region%turbInput%const(MC_SA_RSIG)

! compute vorticities to be used below ----------------------------------------

  CALL TURB_CalcVortic( region )

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv      => region%levels(iLev)%mixt%cv
  tv      => region%levels(iLev)%mixt%tv
  tcv     => region%levels(iLev)%turb%cv
  trhs    => region%levels(iLev)%turb%rhs
  vort    => region%levels(iLev)%turb%vort
  wdist   => region%levels(iLev)%turb%lens
  dsterm  => region%levels(iLev)%turb%dsterm
  vol     => region%levels(iLev)%grid%vol

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
#endif
#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  cv      => region%mixt%cv
  tv      => region%mixt%tv
  tcv     => region%turb%cv
  trhs    => region%turb%rhs
  vort    => region%turb%vort
  wdist   => region%turb%lens
  dsterm  => region%turb%dsterm
  vol     => region%grid%vol

  DO ijkC = ibc,iec
#endif
        vortMag = SQRT( vort(XCOORD,ijkC)*vort(XCOORD,ijkC) + &
                        vort(YCOORD,ijkC)*vort(YCOORD,ijkC) + &
                        vort(ZCOORD,ijkC)*vort(ZCOORD,ijkC) )

        rnuet  = tcv(CV_SA_NUTIL,ijkC)
        nuet   = rnuet/cv(CV_MIXT_DENS,ijkC)
        xi     = rnuet/tv(TV_MIXT_MUEL,ijkC)
        fv1    = xi**3/(xi**3 + cv1**3)
        fv2    = 1._RFREAL - xi/(1._RFREAL + fv1*xi)
        rWdist = 1._RFREAL/wdist(ijkC)
        sTilde = vortMag + nuet*rKappa*rKappa*rWdist*rWdist*fv2
        sTilde = sTilde + REAL_SMALL 
        ro     = nuet*rKappa*rKappa*rWdist*rWdist/sTilde
!        ro     = MIN( ro, 10._RFREAL )
        go     = ro + cw2*(ro**6 - ro)
        fw     = go*((1._RFREAL+cw3**6)/(go**6+cw3**6))**one6th

        prod   = cb1*sTilde*rnuet
        destr  = cv(CV_MIXT_DENS,ijkC)*cw1*fw*(nuet*rWdist)**2
        sTerm  = destr - prod

        trhs(CV_SA_NUTIL,ijkC)   = trhs(CV_SA_NUTIL,ijkC) + vol(ijkC)*sTerm
        dsterm(CV_SA_NUTIL,ijkC) = cb1*(2._RFREAL*sTilde-vortMag) - &
                                   2._RFREAL*cw1*fw*nuet*rWdist*rWdist
#ifdef RFLO
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k
#endif
#ifdef RFLU
  ENDDO       ! ijkC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansSASourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansSASourceTerms.F90,v $
! Revision 1.5  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.8  2004/01/22 03:57:53  wasistho
! fixed dsterm
!
! Revision 1.7  2003/10/25 22:05:43  wasistho
! deactivated (outcommented) clipping on ro
!
! Revision 1.6  2003/10/21 20:31:13  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.5  2003/10/20 20:28:08  wasistho
! made consistent with compressible SA formulation
!
! Revision 1.4  2003/10/16 20:17:38  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
! Revision 1.3  2003/10/14 21:24:54  wasistho
! get minimum of ro and 10.
!
! Revision 1.2  2003/10/07 20:33:46  wasistho
! bug fixed missing nodeOffsets
!
!
!******************************************************************************







