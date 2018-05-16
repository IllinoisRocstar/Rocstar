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
! Purpose: compute diffusion approximation (Rosseland) radiative fluxes by
!          average of variables.
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
! $Id: RADI_ExtinctionCoef.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_ExtinctionCoef( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
#ifdef PLAG
  USE PLAG_ModParameters
#endif
#ifdef PEUL
  USE PEUL_ModParameters
#endif
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, n

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC0, ibc, iec
  INTEGER :: n1, n2, n3, n4, n5, n6, nAddGas, nAddPlag, nAddPeul

  REAL(RFREAL) :: rati, maxVfrac, rTresh
  REAL(RFREAL), POINTER :: coef(:,:), optc(:,:), cellVol(:), goFact(:)
  REAL(RFREAL), POINTER :: vFrac(:,:), pDiam(:,:)
  REAL(RFREAL), ALLOCATABLE :: volFrac(:), extEff(:), rDiam(:)
#ifdef PLAG
  INTEGER :: iPcls, nPcls, nCont, nPclCel
  REAL(RFREAL) :: diaSumCel, volSumCel, volSumPcl
  INTEGER,      POINTER :: pAiv(:,:), pDvPlagVolu(:)
  REAL(RFREAL), POINTER :: pDv(:,:)
#endif
#ifdef PEUL
  INTEGER :: iPtype
  REAL(RFREAL) :: ptDiam, ptEffDens, smokDens
  REAL(RFREAL), POINTER :: peulCv(:,:)
#endif

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_ExtinctionCoef',&
  'RADI_ExtinctionCoef.F90' )

! get dimensions, pointers and constants --------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  coef   => region%levels(iLev)%radi%radCoef
  goFact => region%levels(iLev)%radi%goFact
  optc   => region%radiInput%optConst

  rati = 1.5_RFREAL

! check phase indexing

  IF ((RADI_PHASE_DISPART - RADI_PHASE_GAS /= 1 ) .OR. &
      (RADI_PHASE_CONPART - RADI_PHASE_DISPART /= 1 )) THEN
    CALL ErrorStop( global,ERR_RADI_FIXPARAM,__LINE__, &
         'radiation phase indexing is not consistent' )
  ENDIF

! clear up radiation coefficients before accumulating optical constants

  coef = 0._RFREAL

! extinction coefficient from contribution of AlOxide smoke and Al-dropplets

  IF (region%radiInput%media == RADI_MEDIA_ARTIF) THEN

    ALLOCATE( volFrac(NPHASE), extEff(NPHASE), rDiam(NPHASE) )

    DO n = 1,NPHASE
      volFrac(n) = optc(PHASE_PROP_V,n)
      extEff(n)  = optc(PHASE_PROP_Q,n)
      rDiam(n)   = 1._RFREAL/optc(PHASE_PROP_D,n)
    ENDDO

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

          DO n = 1,NPHASE
            coef(ijkC0,RADI_COEFF_EXTINCT)= coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                                            volFrac(n)*extEff(n)*rDiam(n)
          ENDDO

          coef(ijkC0,RADI_COEFF_EXTINCT) = rati*coef(ijkC0,RADI_COEFF_EXTINCT)
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

    DEALLOCATE( volFrac, extEff, rDiam )

  ELSEIF (region%radiInput%media == RADI_MEDIA_REAL) THEN

    cellVol  => region%levels(iLev)%grid%vol
    nAddGas  = 0
    nAddPlag = 0
    nAddPeul = 0

#ifdef PLAG
    IF (global%plagUsed) THEN
      nAddPlag = 0
    ENDIF
#endif
#ifdef PEUL
    IF (global%peulUsed) THEN
      nAddPeul = region%peulInput%nPtypes - 1
    ENDIF
#endif

    n1 = RADI_PHASE_GAS
    n2 = RADI_PHASE_GAS     + nAddGas
    n3 = RADI_PHASE_DISPART + nAddGas
    n4 = RADI_PHASE_DISPART + nAddGas + nAddPlag
    n5 = RADI_PHASE_CONPART + nAddGas + nAddPlag
    n6 = RADI_PHASE_CONPART + nAddGas + nAddPlag + nAddPeul

    ALLOCATE( vFrac(n6,ibc:iec) )
    ALLOCATE( pDiam(n6,ibc:iec) )

! - initialiaze volume fractions and particles diameters for safety
    vFrac = 0._RFREAL
    pDiam = 1._RFREAL
    vFrac(n1:n2,:) = optc(PHASE_PROP_V,RADI_PHASE_GAS)
    pDiam(n1:n2,:) = optc(PHASE_PROP_D,RADI_PHASE_GAS)

! - disPart volume fraction and diameter
#ifdef PLAG
    IF (global%plagUsed) THEN

      pAiv        => region%levels(iLev)%plag%aiv
      pDv         => region%levels(iLev)%plag%dv
      pDvPlagVolu => region%levels(iLev)%plag%dvPlagVolu

      nPcls    = region%levels(iLev)%plag%nPcls
      nCont    = region%plagInput%nCont

! --- active for non-zero number of particles in region
      IF ( nPcls > 0 ) THEN

        DO k=kpcbeg,kpcend
          DO j=jpcbeg,jpcend
            DO i=ipcbeg,ipcend
              ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

              diaSumCel = 0._RFREAL
              volSumCel = 0._RFREAL
              nPclCel   = 0

              DO iPcls = 1, nPcls

! ------------- only particles in ijkC0 cell
                IF (pAiv(AIV_PLAG_ICELLS,iPcls) == ijkC0) THEN

! --------------- particle volume
                  volSumPcl = SUM(pDv(pDvPlagVolu(:),iPcls))

! --------------- volume of all pcls in this cell
                  volSumCel = volSumCel + volSumPcl

! --------------- sum all diameters and number of particles in this cell
                  diaSumCel = diaSumCel + pDv(DV_PLAG_DIAM,iPcls)
                  nPclCel   = nPclCel + 1
                ENDIF

              ENDDO ! iPcls

              vFrac(n3:n4,ijkC0) = volSumCel/cellVol(ijkC0)
              pDiam(n3:n4,ijkC0) = diaSumCel/REAL(nPclCel,KIND=RFREAL)
            ENDDO ! i
          ENDDO   ! j
        ENDDO     ! k

! ----- copy vFrac and pDiam to dummy cells
        CALL Copy2DumCells( n3,n4,vFrac )
        CALL Copy2DumCells( n3,n4,pDiam )

      ENDIF ! nPcls >0
    ENDIF   ! plagUsed
#endif

! - conPart volume fraction and diameter
#ifdef PEUL
    IF (global%peulUsed) THEN

      peulCv => region%levels(iLev)%peul%cv

! --- there are nPtypes smoke particles, each treated as different phase
      DO iPtype = 1,region%peulInput%nPtypes

        ptDiam    = region%peulInput%ptypes(iPtype)%diam
        ptEffDens = region%peulInput%ptypes(iPtype)%denseff

        DO k=kdcbeg,kdcend
          DO j=jdcbeg,jdcend
            DO i=idcbeg,idcend
              ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

              smokDens = MAX( peulCv(CV_PEUL_DENS+iPtype-1,ijkC0), 0._RFREAL )
              vFrac(n5+iPtype-1,ijkC0) = smokDens/ptEffDens
              pDiam(n5+iPtype-1,ijkC0) = ptDiam

            ENDDO ! i
          ENDDO   ! j
        ENDDO     ! k
      ENDDO       ! iPtype
    ENDIF         ! peulUsed
#endif

    ALLOCATE( volFrac(n6) )
    ALLOCATE( extEff(n6) )
    volFrac = 0._RFREAL
    extEff  = 0._RFREAL

    extEff(n1:n2)  = optc(PHASE_PROP_Q,RADI_PHASE_GAS)
    extEff(n3:n4)  = optc(PHASE_PROP_Q,RADI_PHASE_DISPART)
    extEff(n5:n6)  = optc(PHASE_PROP_Q,RADI_PHASE_CONPART)

    maxVfrac = 0._RFREAL

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)

          DO n = n1, n2
            volFrac(n) = MIN( optc(PHASE_PROP_V,RADI_PHASE_GAS), &
                              vFrac(n,ijkC0) )
          ENDDO
          DO n = n3, n4
            volFrac(n) = MIN( optc(PHASE_PROP_V,RADI_PHASE_DISPART), &
                              vFrac(n,ijkC0) )
          ENDDO
          DO n = n5, n6
            volFrac(n) = MIN( optc(PHASE_PROP_V,RADI_PHASE_CONPART), &
                              vFrac(n,ijkC0) )
          ENDDO

          maxVfrac = MAX( maxVfrac, SUM( volFrac(:) ) )

          DO n = n1,n6
            coef(ijkC0,RADI_COEFF_EXTINCT)= coef(ijkC0,RADI_COEFF_EXTINCT)+ &
                                            volFrac(n)*extEff(n)/pDiam(n,ijkC0)
          ENDDO

          coef(ijkC0,RADI_COEFF_EXTINCT) = rati*coef(ijkC0,RADI_COEFF_EXTINCT)
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

    DEALLOCATE( vFrac, pDiam, volFrac, extEff )
    IF (maxVfrac > 1._RFREAL) GOTO 10

  ENDIF         ! media real

! compute go-ahead factor based on EC treshold to obey diffusion approximation

  rTresh = 1._RFREAL/RADI_REAL_ECMIN

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC0 = IndIJK(i ,j ,k ,iCOff,ijCOff)
        goFact(ijkC0) = AINT( coef(ijkC0,RADI_COEFF_EXTINCT)*rTresh )
        goFact(ijkC0) = MIN( goFact(ijkC0), 1._RFREAL )
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_RADI_MULPHASE,__LINE__, &
       'total volume fractions of particles and smoke > 1.' )

999 CONTINUE

  CALL DeregisterFunction( region%global )

! =============================================================================
!   Copy given vector variables to dummy cells
! =============================================================================

CONTAINS

  SUBROUTINE Copy2DumCells( m1,m2,var )

! ... parameters
    INTEGER :: m1, m2, ibegc, iendc
    REAL(RFREAL), POINTER :: var(:,:)

! ... local variables
    INTEGER :: l, ijkC, ijkCi

! - I direction

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend

        DO i=idcbeg,ipcbeg-1
          ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkCi = IndIJK(1     ,j     ,k     ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO
        DO i=ipcend+1,idcend
          ijkC  = IndIJK(i      ,j     ,k     ,iCOff,ijCOff)
          ijkCi = IndIJK(ipcend ,j     ,k     ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

! - J direction

    DO k=kdcbeg,kdcend
      DO i=idcbeg,idcend

        DO j=jdcbeg,jpcbeg-1
          ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkCi = IndIJK(i     ,1     ,k     ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO

        DO j=jpcend+1,jdcend
          ijkC  = IndIJK(i     ,j      ,k     ,iCOff,ijCOff)
          ijkCi = IndIJK(i     ,jpcend ,k     ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO   ! j
      ENDDO     ! i
    ENDDO       ! k

! - K direction

    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend

        DO k=kdcbeg,kpcbeg-1
          ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkCi = IndIJK(i     ,j     ,1     ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO

        DO k=kpcend+1,kdcend
          ijkC  = IndIJK(i     ,j     ,k      ,iCOff,ijCOff)
          ijkCi = IndIJK(i     ,j     ,kpcend ,iCOff,ijCOff)

          DO l = m1, m2
            var(l,ijkC) = var(l,ijkCi)
          ENDDO
        ENDDO   ! k
      ENDDO     ! i
    ENDDO       ! j

  END SUBROUTINE Copy2DumCells

END SUBROUTINE RADI_ExtinctionCoef

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_ExtinctionCoef.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.3  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2003/07/30 22:23:02  wasistho
! enter part and smoke data into radiation
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!******************************************************************************







