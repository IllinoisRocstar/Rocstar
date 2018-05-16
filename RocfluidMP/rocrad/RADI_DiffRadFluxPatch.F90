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
!          through a patch by using an average of variables.
!
! Description: none.
!
! Input: region = data of current region
!        patch  = current patch.
!
! Output: region%levels%mixt%rhs = radiative fluxes added to the residual.
!
! Notes: In addition to the energy source term due to radiation (div(qr)), 
!        radiant fluxes normal to surface (cell faces) are stored in data 
!        structure to be delivered later to combustion module, through Genx,
!        for burn-rate calculations. The radiant fluxes are stored always in
!        positive directions w.r.t lbound = 1, 3, and 5. Opposite sign 
!        should be applied later for lbound = 2, 4, and 6 before transfering
!        the radiant fluxes through Genx.
!
!******************************************************************************
!
! $Id: RADI_DiffRadFluxPatch.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_DiffRadFluxPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, n1, n2

! ... local variables
  INTEGER :: iLev, lbound, bcType, distrib, flowModel, bcOpt
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCD, ijkCB0, ijkCB1, ijkNB
  INTEGER :: inode, jnode, knode
  INTEGER :: idbeg, idend, jdbeg, jdend, kdbeg, kdend, i2d, nOff

  REAL(RFREAL)          :: sgn, gFaca, flima, coefa, tempa, tBurn, tWall, mRate
  REAL(RFREAL)          :: stBoltz, rati, beta, qrx, qry, qrz, fr, modSf, sf(3)
  REAL(RFREAL), POINTER :: rhs(:,:), sFace(:,:), vals(:,:), grad(:,:) 
  REAL(RFREAL), POINTER :: rdv(:,:), coef(:,:), goFact(:), flim(:), qr(:)
  REAL(RFREAL), POINTER :: wvInt(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_DiffRadFluxPatch',&
  'RADI_DiffRadFluxPatch.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  bcType    = patch%bcType
  nOff      = ABS(patch%l1end-patch%l1beg) + 1
  distrib   = patch%mixt%distrib
  flowModel = region%mixtInput%flowModel

  rhs   => region%levels(iLev)%mixt%rhs
  vals  => patch%mixt%vals
  rdv   => region%levels(iLev)%radi%dv
  flim  => region%levels(iLev)%radi%fluxLim
  coef  => region%levels(iLev)%radi%radCoef
  goFact=> region%levels(iLev)%radi%goFact
  wvInt => region%levels(iLev)%radi%wvInt

! take the right face vector and make it point outwards

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

! get the appropriate face vector, gradients and radiant flux

  IF (lbound==1 .OR. lbound==2) THEN
    sFace => region%levels(iLev)%grid%si
    grad  => region%levels(iLev)%mixt%gradi
    qr    => region%levels(iLev)%radi%qri
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    sFace => region%levels(iLev)%grid%sj
    grad  => region%levels(iLev)%mixt%gradj
    qr    => region%levels(iLev)%radi%qrj
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    sFace => region%levels(iLev)%grid%sk
    grad  => region%levels(iLev)%mixt%gradk
    qr    => region%levels(iLev)%radi%qrk
  ENDIF
 
! get radiation and rk-stage constants

  stBoltz = region%radiInput%stBoltz
  rati    = -16._RFREAL
  beta    = region%mixtInput%betrk(region%irkStep)

! stationary grid -------------------------------------------------------------
! slip wall

  IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) THEN

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCB1 = IndIJK(i+idir ,j+jdir ,k+kdir ,iCOff,ijCOff)  ! interior
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

          gFaca  = goFact(ijkCB0) 
          flima  = 0.5_RFREAL*(3._RFREAL*flim(ijkCB0)-flim(ijkCB1))
          tempa  = 0.5_RFREAL*(3._RFREAL*rdv(DV_RADI_TEFF,ijkCB0)- &
                                         rdv(DV_RADI_TEFF,ijkCB1))
          coefa  = 0.5_RFREAL*(3._RFREAL*coef(ijkCB0,RADI_COEFF_EXTINCT)- &
                                         coef(ijkCB1,RADI_COEFF_EXTINCT))

          coefa = gFaca*rati*flima*stBoltz*tempa**3/coefa

          qrx = coefa*grad(GR_MIXT_TX,ijkNB)
          qry = coefa*grad(GR_MIXT_TY,ijkNB)
          qrz = coefa*grad(GR_MIXT_TZ,ijkNB)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

          fr = qrx*sf(1) + qry*sf(2) + qrz*sf(3)

          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) - fr*beta

! ------- store rad. flux normal to surf. in positive direction (hence sgn* )
          modSf     = SQRT( sf(1)*sf(1) + sf(2)*sf(2) + sf(3)*sf(3) )  
          qr(ijkNB) = -sgn*fr/modSf
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! noslip wall

  ELSE IF (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

    bcOpt = patch%mixt%switches(BCSWI_NOSLIP_ADIABAT)

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

! ------- adiabatic wall
          IF (bcOpt == BCOPT_ADIABAT) THEN
            tempa  = rdv(DV_RADI_TEFF,ijkCB0)
          ELSE

! ------- prescribed wall temperature
            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d   = distrib * IndIJ(n1,n2,nOff)
            tWall = vals(BCDAT_NOSLIP_TWALL,i2d)
            tempa = tWall
          ENDIF

          gFaca  = goFact(ijkCB0) 
          flima  = flim(ijkCB0) 
          coefa  = coef(ijkCB0,RADI_COEFF_EXTINCT)
          coefa  = gFaca*rati*flima*stBoltz*tempa**3/coefa

          qrx = coefa*grad(GR_MIXT_TX,ijkNB)
          qry = coefa*grad(GR_MIXT_TY,ijkNB)
          qrz = coefa*grad(GR_MIXT_TZ,ijkNB)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

          fr = qrx*sf(1) + qry*sf(2) + qrz*sf(3)

          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) - fr*beta

! ------- store rad. flux normal to surf. in positive direction (hence sgn* )
          modSf     = SQRT( sf(1)*sf(1) + sf(2)*sf(2) + sf(3)*sf(3) )  
          qr(ijkNB) = -sgn*fr/modSf
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! injection boundary 

  ELSE IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
           

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCB1 = IndIJK(i+idir ,j+jdir ,k+kdir ,iCOff,ijCOff)  ! interior
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

          IF      (lbound==1 .OR. lbound==2) THEN
            n1 = j - jbeg
            n2 = k - kbeg
          ELSE IF (lbound==3 .OR. lbound==4) THEN
            n1 = k - kbeg
            n2 = i - ibeg
          ELSE IF (lbound==5 .OR. lbound==6) THEN
            n1 = i - ibeg
            n2 = j - jbeg
          ENDIF
          i2d    = distrib * IndIJ(n1,n2,nOff)
          mRate  = vals(BCDAT_INJECT_MFRATE,i2d)
          tBurn  = vals(BCDAT_INJECT_TEMP  ,i2d)
          gFaca  = goFact(ijkCB0)

          IF (mRate > 0._RFREAL) THEN        ! surface burning
            flima = flim(ijkCB0)
            tempa = tBurn
            coefa = coef(ijkCB0,RADI_COEFF_EXTINCT)

          ELSE                               ! not burning - slip/noslip wall
            IF (flowModel == FLOW_EULER) THEN
              flima = 0.5_RFREAL*(3._RFREAL*flim(ijkCB0)-flim(ijkCB1))
              tempa = 0.5_RFREAL*(3._RFREAL*rdv(DV_RADI_TEFF,ijkCB0)- &
                                            rdv(DV_RADI_TEFF,ijkCB1))
              coefa = 0.5_RFREAL*(3._RFREAL*coef(ijkCB0,RADI_COEFF_EXTINCT)- &
                                            coef(ijkCB1,RADI_COEFF_EXTINCT))

            ELSE                             ! treated as NS adiabatic
              flima = flim(ijkCB0)
              tempa = rdv(DV_RADI_TEFF,ijkCB0)
              coefa = coef(ijkCB0,RADI_COEFF_EXTINCT)
            ENDIF
          ENDIF

          coefa  = gFaca*rati*flima*stBoltz*tempa**3/coefa

          qrx = coefa*grad(GR_MIXT_TX,ijkNB)
          qry = coefa*grad(GR_MIXT_TY,ijkNB)
          qrz = coefa*grad(GR_MIXT_TZ,ijkNB)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

          fr = qrx*sf(1) + qry*sf(2) + qrz*sf(3)

          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) - fr*beta

! ------- store rad. flux normal to surf. in positive direction (hence sgn* )
          modSf     = SQRT( sf(1)*sf(1) + sf(2)*sf(2) + sf(3)*sf(3) )  
          qr(ijkNB) = -sgn*fr/modSf
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! non-conforming region interface

  ELSE IF (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! everything else

  ELSE

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
          ijkCD  = IndIJK(i-idir ,j-jdir ,k-kdir ,iCOff,ijCOff)  ! dummy
          ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

          gFaca  = 0.5_RFREAL*(goFact(ijkCB0)+goFact(ijkCD)) 
          flima  = 0.5_RFREAL*(flim(ijkCB0)+flim(ijkCD)) 
          tempa  = 0.5_RFREAL*(rdv(DV_RADI_TEFF,ijkCB0) + &
                               rdv(DV_RADI_TEFF,ijkCD))
          coefa  = 0.5_RFREAL*(coef(ijkCB0,RADI_COEFF_EXTINCT) + &
                               coef(ijkCD ,RADI_COEFF_EXTINCT))

          coefa = gFaca*rati*flima*stBoltz*tempa**3/coefa

          qrx = coefa*grad(GR_MIXT_TX,ijkNB)
          qry = coefa*grad(GR_MIXT_TY,ijkNB)
          qrz = coefa*grad(GR_MIXT_TZ,ijkNB)

          sf(1)  = sgn*sFace(XCOORD,ijkNB)
          sf(2)  = sgn*sFace(YCOORD,ijkNB)
          sf(3)  = sgn*sFace(ZCOORD,ijkNB)

          fr = qrx*sf(1) + qry*sf(2) + qrz*sf(3)

          rhs(CV_MIXT_ENER,ijkCB0) = rhs(CV_MIXT_ENER,ijkCB0) - fr*beta

! ------- store rad. flux normal to surf. in positive direction (hence sgn* )
          modSf     = SQRT( sf(1)*sf(1) + sf(2)*sf(2) + sf(3)*sf(3) )  
          qr(ijkNB) = -sgn*fr/modSf
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  ENDIF         ! bcType

! add contribution of grad(T) from patch for intensity calculation later

  IF (lbound==1 .OR. lbound==2) THEN
    idbeg = ibeg - idir*region%nDumCells
    idend = iend - idir*region%nDumCells
    jdbeg = jbeg - region%nDumCells
    jdend = jend + region%nDumCells
    kdbeg = kbeg - region%nDumCells
    kdend = kend + region%nDumCells

  ELSE IF (lbound==3 .OR. lbound==4) THEN
    jdbeg = jbeg - jdir*region%nDumCells
    jdend = jend - jdir*region%nDumCells
    idbeg = ibeg - region%nDumCells
    idend = iend + region%nDumCells
    kdbeg = kbeg - region%nDumCells
    kdend = kend + region%nDumCells

  ELSE IF (lbound==5 .OR. lbound==6) THEN
    kdbeg = kbeg - kdir*region%nDumCells
    kdend = kend - kdir*region%nDumCells
    jdbeg = jbeg - region%nDumCells
    jdend = jend + region%nDumCells
    idbeg = ibeg - region%nDumCells
    idend = iend + region%nDumCells
  ENDIF

  DO k=kdbeg,kdend
    DO j=jdbeg,jdend
      DO i=idbeg,idend
        ijkCB0 = IndIJK(i      ,j      ,k      ,iCOff,ijCOff)  ! boundary
        ijkNB  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

        wvInt(XCOORD,ijkCB0) = wvInt(XCOORD,ijkCB0) + grad(GR_MIXT_TX,ijkNB)
        wvInt(YCOORD,ijkCB0) = wvInt(YCOORD,ijkCB0) + grad(GR_MIXT_TY,ijkNB)
        wvInt(ZCOORD,ijkCB0) = wvInt(ZCOORD,ijkCB0) + grad(GR_MIXT_TZ,ijkNB)
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_DiffRadFluxPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_DiffRadFluxPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:11  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.5  2004/09/22 01:32:10  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.4  2004/09/18 17:41:54  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.3  2003/08/13 21:18:12  wasistho
! corrected sign and rk-stage of radiation flux
!
! Revision 1.2  2003/07/30 22:24:34  wasistho
! enter part and smoke data into radiation
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!******************************************************************************







