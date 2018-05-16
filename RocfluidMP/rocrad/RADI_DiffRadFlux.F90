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
! Purpose: compute diffusion approximation, FLDSRC or ROSS, radiative fluxes 
!          by average of variables.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = radiative fluxes added to the residual.
!
! Notes: grad(T) is averaged to cell center in wvInt for computation of
!        radiation intensity later.
!
!******************************************************************************
!
! $Id: RADI_DiffRadFlux.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_DiffRadFlux( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE RADI_ModInterfaces, ONLY : RADI_DiffRadFluxPatch
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iPatch

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC0, ijkC1, ijkN

  REAL(RFREAL)          :: gFaca, flima, coefa, tempa, stBoltz, rati, beta
  REAL(RFREAL)          :: qrx, qry, qrz, modSf, fr
  REAL(RFREAL), POINTER :: rhs(:,:), si(:,:), sj(:,:), sk(:,:)
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  REAL(RFREAL), POINTER :: rdv(:,:), qri(:), qrj(:), qrk(:) 
  REAL(RFREAL), POINTER :: goFact(:), flim(:), coef(:,:), wvInt(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_DiffRadFlux',&
  'RADI_DiffRadFlux.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk

  rhs    => region%levels(iLev)%mixt%rhs
  gradi  => region%levels(iLev)%mixt%gradi
  gradj  => region%levels(iLev)%mixt%gradj
  gradk  => region%levels(iLev)%mixt%gradk

  rdv    => region%levels(iLev)%radi%dv
  flim   => region%levels(iLev)%radi%fluxLim
  coef   => region%levels(iLev)%radi%radCoef
  goFact => region%levels(iLev)%radi%goFact
  qri    => region%levels(iLev)%radi%qri
  qrj    => region%levels(iLev)%radi%qrj
  qrk    => region%levels(iLev)%radi%qrk
  wvInt  => region%levels(iLev)%radi%wvInt
 
  wvInt   = 0._RFREAL
  stBoltz = region%radiInput%stBoltz
  rati    = -16._RFREAL
  beta    = region%mixtInput%betrk(region%irkStep)

! flux in i-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg+1,ipcend
        ijkC0 = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)

        gFaca = 0.5_RFREAL*(goFact(ijkC0)+goFact(ijkC1)) 
        flima = 0.5_RFREAL*(flim(ijkC0)+flim(ijkC1)) 
        coefa = 0.5_RFREAL*(coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                            coef(ijkC1,RADI_COEFF_EXTINCT))
        tempa = 0.5_RFREAL*(rdv(DV_RADI_TEFF,ijkC0)+rdv(DV_RADI_TEFF,ijkC1))
        coefa = gFaca*rati*flima*stBoltz*tempa**3/coefa

        qrx   = coefa*gradi(GR_MIXT_TX,ijkN)
        qry   = coefa*gradi(GR_MIXT_TY,ijkN)
        qrz   = coefa*gradi(GR_MIXT_TZ,ijkN)

        fr    = qrx*si(XCOORD,ijkN)+ &
                qry*si(YCOORD,ijkN)+ &
                qrz*si(ZCOORD,ijkN)

        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) - fr*beta
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) + fr*beta

! ----- store radiant flux normal to cell face in positive direction in qri
        modSf     = SQRT( si(XCOORD,ijkN)*si(XCOORD,ijkN) + &
                          si(YCOORD,ijkN)*si(YCOORD,ijkN) + &
                          si(ZCOORD,ijkN)*si(ZCOORD,ijkN))
        qri(ijkN) = -fr/modSf
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg+1,idcend
        ijkC0 = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkN  = IndIJK(i  ,j,k,iNOff,ijNOff)

        wvInt(XCOORD,ijkC0) = wvInt(XCOORD,ijkC0) + gradi(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC0) = wvInt(YCOORD,ijkC0) + gradi(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC0) = wvInt(ZCOORD,ijkC0) + gradi(GR_MIXT_TZ,ijkN)

        wvInt(XCOORD,ijkC1) = wvInt(XCOORD,ijkC1) + gradi(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC1) = wvInt(YCOORD,ijkC1) + gradi(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC1) = wvInt(ZCOORD,ijkC1) + gradi(GR_MIXT_TZ,ijkN)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! flux in j-direction (except through boundary) -------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg+1,jpcend
      DO i=ipcbeg,ipcend
        ijkC0 = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)

        gFaca = 0.5_RFREAL*(goFact(ijkC0)+goFact(ijkC1)) 
        flima = 0.5_RFREAL*(flim(ijkC0)+flim(ijkC1)) 
        coefa = 0.5_RFREAL*(coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                            coef(ijkC1,RADI_COEFF_EXTINCT))
        tempa = 0.5_RFREAL*(rdv(DV_RADI_TEFF,ijkC0)+rdv(DV_RADI_TEFF,ijkC1))

        coefa = gFaca*rati*flima*stBoltz*tempa**3/coefa

        qrx   = coefa*gradj(GR_MIXT_TX,ijkN)
        qry   = coefa*gradj(GR_MIXT_TY,ijkN)
        qrz   = coefa*gradj(GR_MIXT_TZ,ijkN)

        fr    = qrx*sj(XCOORD,ijkN)+ &
                qry*sj(YCOORD,ijkN)+ &
                qrz*sj(ZCOORD,ijkN)

        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) - fr*beta
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) + fr*beta

! ----- store radiant flux normal to cell face in positive direction in qrj
        modSf     = SQRT( sj(XCOORD,ijkN)*sj(XCOORD,ijkN) + &
                          sj(YCOORD,ijkN)*sj(YCOORD,ijkN) + &
                          sj(ZCOORD,ijkN)*sj(ZCOORD,ijkN))
        qrj(ijkN) = -fr/modSf
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO k=kdcbeg,kdcend
    DO j=jdcbeg+1,jdcend
      DO i=idcbeg,idcend
        ijkC0 = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkN  = IndIJK(i,j  ,k,iNOff,ijNOff)

        wvInt(XCOORD,ijkC0) = wvInt(XCOORD,ijkC0) + gradj(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC0) = wvInt(YCOORD,ijkC0) + gradj(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC0) = wvInt(ZCOORD,ijkC0) + gradj(GR_MIXT_TZ,ijkN)

        wvInt(XCOORD,ijkC1) = wvInt(XCOORD,ijkC1) + gradj(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC1) = wvInt(YCOORD,ijkC1) + gradj(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC1) = wvInt(ZCOORD,ijkC1) + gradj(GR_MIXT_TZ,ijkN)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! flux in k-direction (except through boundary) -------------------------------

  DO k=kpcbeg+1,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkC0 = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)

        gFaca = 0.5_RFREAL*(goFact(ijkC0)+goFact(ijkC1)) 
        flima = 0.5_RFREAL*(flim(ijkC0)+flim(ijkC1)) 
        coefa = 0.5_RFREAL*(coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                            coef(ijkC1,RADI_COEFF_EXTINCT))
        tempa = 0.5_RFREAL*(rdv(DV_RADI_TEFF,ijkC0)+rdv(DV_RADI_TEFF,ijkC1))

        coefa = gFaca*rati*flima*stBoltz*tempa**3/coefa

        qrx   = coefa*gradk(GR_MIXT_TX,ijkN)
        qry   = coefa*gradk(GR_MIXT_TY,ijkN)
        qrz   = coefa*gradk(GR_MIXT_TZ,ijkN)

        fr    = qrx*sk(XCOORD,ijkN)+ &
                qry*sk(YCOORD,ijkN)+&
                qrz*sk(ZCOORD,ijkN)

        rhs(CV_MIXT_ENER,ijkC0) = rhs(CV_MIXT_ENER,ijkC0) - fr*beta
        rhs(CV_MIXT_ENER,ijkC1) = rhs(CV_MIXT_ENER,ijkC1) + fr*beta

! ----- store radiant flux normal to cell face in positive direction in qrk
        modSf     = SQRT( sk(XCOORD,ijkN)*sk(XCOORD,ijkN) + &
                          sk(YCOORD,ijkN)*sk(YCOORD,ijkN) + &
                          sk(ZCOORD,ijkN)*sk(ZCOORD,ijkN))
        qrk(ijkN) = -fr/modSf
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

  DO k=kdcbeg+1,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC0 = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkN  = IndIJK(i,j,k  ,iNOff,ijNOff)

        wvInt(XCOORD,ijkC0) = wvInt(XCOORD,ijkC0) + gradk(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC0) = wvInt(YCOORD,ijkC0) + gradk(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC0) = wvInt(ZCOORD,ijkC0) + gradk(GR_MIXT_TZ,ijkN)

        wvInt(XCOORD,ijkC1) = wvInt(XCOORD,ijkC1) + gradk(GR_MIXT_TX,ijkN)
        wvInt(YCOORD,ijkC1) = wvInt(YCOORD,ijkC1) + gradk(GR_MIXT_TY,ijkN)
        wvInt(ZCOORD,ijkC1) = wvInt(ZCOORD,ijkC1) + gradk(GR_MIXT_TZ,ijkN)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fluxes through boundaries ---------------------------------------------------

  DO iPatch=1,region%nPatches
    CALL RADI_DiffRadFluxPatch( region,region%levels(iLev)%patches(iPatch) )
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_DiffRadFlux

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_DiffRadFlux.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.5  2004/09/22 01:32:03  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.4  2004/09/18 17:41:48  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.3  2003/08/13 21:18:07  wasistho
! corrected sign and rk-stage of radiation flux
!
! Revision 1.2  2003/07/30 22:24:28  wasistho
! enter part and smoke data into radiation
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!******************************************************************************







