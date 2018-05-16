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
! Purpose: compute the radiation flux limiter based on the model selected.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%fluxLim = flux limiter in diff. based radiation
!
! Notes: none. 
!
!******************************************************************************
!
! $Id: RADI_FluxLimiter.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FluxLimiter( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: ijkC, ijkN, ijkNi, ijkNj, ijkNk, errFl

  REAL(RFREAL) :: one3rd, one6th, rparam, gradEMod, coefc, radEnerg
  REAL(RFREAL), POINTER :: dv(:,:), rcv(:,:), rdv(:,:), coef(:,:), fluxLim(:) 
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  REAL(RFREAL), POINTER :: rGradi(:,:), rGradj(:,:), rGradk(:,:)
  REAL(RFREAL), ALLOCATABLE :: gradc(:,:,:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FluxLimiter',&
  'RADI_FluxLimiter.F90' )

! get constants, dimensions and pointers --------------------------------------

  one3rd = 1._RFREAL/3._RFREAL
  one6th = 1._RFREAL/6._RFREAL
  iLev   = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  gradi  => region%levels(iLev)%mixt%gradi
  gradj  => region%levels(iLev)%mixt%gradj
  gradk  => region%levels(iLev)%mixt%gradk
  dv     => region%levels(iLev)%mixt%dv
  rdv    => region%levels(iLev)%radi%dv
  coef   => region%levels(iLev)%radi%radCoef

  IF (region%radiInput%radiModel == RADI_MODEL_FLDSRC .OR. & 
      region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    ALLOCATE( gradc(idcbeg:idcend,jdcbeg:jdcend,kdcbeg:kdcend,3),stat=errFl )
    region%global%error = errFl
    CALL ErrorStop( region%global,ERR_ALLOCATE,__LINE__ )
  ENDIF
 
! get effective temperature --------------------------------------------------

  IF ((region%radiInput%radiModel == RADI_MODEL_ROSS) .OR. &
      (region%radiInput%radiModel == RADI_MODEL_FLDSRC .AND. & 
       region%radiInput%fluxLim == FLD_LIM_NONE)) THEN
    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)
          fluxLim(ijkC) = one3rd
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
  ELSEIF ((region%radiInput%radiModel == RADI_MODEL_FLDSRC) .AND. &
          (region%radiInput%fluxLim == FLD_LIM_LP)) THEN

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkN  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkNi = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNj = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkNk = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          gradc(i,j,k,1) = one6th* &
                         (gradi(GR_MIXT_TX,ijkN)+gradi(GR_MIXT_TX,ijkNi) + &
                          gradj(GR_MIXT_TX,ijkN)+gradj(GR_MIXT_TX,ijkNj) + &
                          gradk(GR_MIXT_TX,ijkN)+gradk(GR_MIXT_TX,ijkNk))
          gradc(i,j,k,2) = one6th* &
                         (gradi(GR_MIXT_TY,ijkN)+gradi(GR_MIXT_TY,ijkNi) + &
                          gradj(GR_MIXT_TY,ijkN)+gradj(GR_MIXT_TY,ijkNj) + &
                          gradk(GR_MIXT_TY,ijkN)+gradk(GR_MIXT_TY,ijkNk))
          gradc(i,j,k,3) = one6th* &
                         (gradi(GR_MIXT_TZ,ijkN)+gradi(GR_MIXT_TZ,ijkNi) + &
                          gradj(GR_MIXT_TZ,ijkN)+gradj(GR_MIXT_TZ,ijkNj) + &
                          gradk(GR_MIXT_TZ,ijkN)+gradk(GR_MIXT_TZ,ijkNk))
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)
          gradEMod = SQRT( gradc(i,j,k,1)*gradc(i,j,k,1) + &
                           gradc(i,j,k,2)*gradc(i,j,k,2) + &
                           gradc(i,j,k,3)*gradc(i,j,k,3) )
          coefc    = coef(ijkC,RADI_COEFF_EXTINCT)
          radEnerg = ABS( dv(DV_MIXT_TEMP,ijkC) )
          rparam   = gradEMod/(coefc*radEnerg)
          fluxLim(ijkC) = fLimiterLP(rparam)
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
  ELSEIF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN

    rcv    => region%levels(iLev)%radi%cv
    rGradi => region%levels(iLev)%radi%gradi
    rGradj => region%levels(iLev)%radi%gradj
    rGradk => region%levels(iLev)%radi%gradk

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkN  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkNi = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNj = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkNk = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          gradc(i,j,k,1) = one6th* &
                         (rGradi(GR_RADI_EX,ijkN)+rGradi(GR_RADI_EX,ijkNi) + &
                          rGradj(GR_RADI_EX,ijkN)+rGradj(GR_RADI_EX,ijkNj) + &
                          rGradk(GR_RADI_EX,ijkN)+rGradk(GR_RADI_EX,ijkNk))
          gradc(i,j,k,2) = one6th* &
                         (rGradi(GR_RADI_EY,ijkN)+rGradi(GR_RADI_EY,ijkNi) + &
                          rGradj(GR_RADI_EY,ijkN)+rGradj(GR_RADI_EY,ijkNj) + &
                          rGradk(GR_RADI_EY,ijkN)+rGradk(GR_RADI_EY,ijkNk))
          gradc(i,j,k,3) = one6th* &
                         (rGradi(GR_RADI_EZ,ijkN)+rGradi(GR_RADI_EZ,ijkNi) + &
                          rGradj(GR_RADI_EZ,ijkN)+rGradj(GR_RADI_EZ,ijkNj) + &
                          rGradk(GR_RADI_EZ,ijkN)+rGradk(GR_RADI_EZ,ijkNk))
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k

    DO k=kdcbeg,kdcend
      DO j=jdcbeg,jdcend
        DO i=idcbeg,idcend
          ijkC = IndIJK(i  ,j,k,iCOff,ijCOff)
          gradEMod = SQRT( gradc(i,j,k,1)*gradc(i,j,k,1) + &
                           gradc(i,j,k,2)*gradc(i,j,k,2) + &
                           gradc(i,j,k,3)*gradc(i,j,k,3) )
          coefc    = coef(ijkC,RADI_COEFF_EXTINCT)
          radEnerg = rcv(CV_RADI_ENER,ijkC)
          rparam   = gradEMod/(coefc*radEnerg)
          fluxLim(ijkC) = fLimiterLP(rparam)
        ENDDO   ! i
      ENDDO   ! j
    ENDDO   ! k
  ENDIF

! deallocate arrays

  IF (ALLOCATED( gradc )) THEN
    DEALLOCATE( gradc, stat=errFl )
    region%global%error = errFl
    CALL ErrorStop( region%global,ERR_DEALLOCATE,__LINE__ )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! #############################################################################

  CONTAINS

    REAL(RFREAL) FUNCTION fLimiterLP( r )

      REAL(RFREAL) :: r

      flimiterLP = (2._RFREAL + r)/(6._RFREAL + 3._RFREAL*r + r*r)

    END FUNCTION fLimiterLP

END SUBROUTINE RADI_FluxLimiter

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FluxLimiter.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.2  2004/09/22 01:32:19  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.1  2004/09/18 18:02:07  wasistho
! initial import for Limited Flux Diffusion radiation
!
!
!
!******************************************************************************







