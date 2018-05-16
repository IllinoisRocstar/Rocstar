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
! Purpose: Write out RADI user input and parameters derived from it for
!          checking purposes.
!
! Description: none.
!
! Input: region = user input in current region.
!
! Output: to standard output (monitor).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_PrintUserInput.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_PrintUserInput( region )
#endif
#ifdef RFLU
SUBROUTINE RADI_PrintUserInput
#endif

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModRadiation
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region) :: region
#endif

! ... loop variables
  INTEGER :: n

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: radiModel, media, solMethod, nPol, nAzi, nAng
  INTEGER :: fluxLim, discrType, discrOrder
  REAL(RFREAL) :: extCoef, smoocf, k2, ook4
  REAL(RFREAL), POINTER :: angles(:,:), optConst(:,:)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_PrintUserInput',&
  'RADI_PrintUserInput.F90' )

! input parameters ------------------------------------------------------------

#ifdef RFLO
  radiModel =  region%radiInput%radiModel
  media     =  region%radiInput%media
  fluxLim   =  region%radiInput%fluxLim
  smoocf    =  region%radiInput%smoocf
  discrType =  region%radiInput%spaceDiscr
  k2        =  region%radiInput%vis2
  ook4      =  1._RFREAL/region%radiInput%vis4
  discrOrder=  region%radiInput%spaceOrder
  solMethod =  region%radiInput%solMethod
  nPol      =  region%radiInput%nPol
  nAzi      =  region%radiInput%nAzi
  nAng      =  region%radiInput%nAng
  angles    => region%radiInput%angles
  optConst  => region%radiInput%optConst
#endif
#ifdef RFLU
  radiModel =  radiInput%radiModel
  media     =  radiInput%media
  fluxLim   =  radiInput%fluxLim
  smoocf    =  radiInput%smoocf
  discrType =  radiInput%spaceDiscr
  k2        =  radiInput%vis2
  ook4      =  1._RFREAL/radiInput%vis4
  discrOrder=  radiInput%spaceOrder
  solMethod =  radiInput%solMethod
  nPol      =  radiInput%nPol
  nAzi      =  radiInput%nAzi
  nAng      =  radiInput%nAng
  angles    => radiInput%angles
  optConst  => radiInput%optConst
#endif

! output selected general parameters to monitor

  WRITE(STDOUT,1005) SOLVER_NAME//' Radiation:'

  IF (radiModel == RADI_MODEL_ROSS) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation model = Rosseland Diff. approx.'
  ELSEIF (radiModel == RADI_MODEL_FLDSRC) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation model = simple Limited Flux Diff.'
  ELSEIF (radiModel == RADI_MODEL_FLDTRAN) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation model = full Limited Flux Diff.'
  ELSEIF (radiModel == RADI_MODEL_RTEGRAY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation model = RTE gray gas'
  ELSEIF (radiModel == RADI_MODEL_RTEBAND) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation model = RTE non-gray gas'
  ENDIF

  IF (media == RADI_MEDIA_ARTIF) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiative media = artificial multi-phase'
  ELSEIF (media == RADI_MEDIA_REAL) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' radiative media = real multi-phase'
  ENDIF

  WRITE(STDOUT,1030) SOLVER_NAME//' relevant input-data based on above model:'

! optical constants and extinction coefficients

  IF (media == RADI_MEDIA_ARTIF) THEN
    IF (optConst(PHASE_PROP_V,RADI_PHASE_GAS) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase     = gas'
      WRITE(STDOUT,1010) SOLVER_NAME//'   vol.fraction  =', &
                         optConst(PHASE_PROP_V,RADI_PHASE_GAS)
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction  =', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_GAS)
      WRITE(STDOUT,1010) SOLVER_NAME//'   molecule diam.=', &
                         optConst(PHASE_PROP_D,RADI_PHASE_GAS)
    ENDIF
    IF (optConst(PHASE_PROP_V,RADI_PHASE_CONPART) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase     = continuum particles'
      WRITE(STDOUT,1010) SOLVER_NAME//'   vol.fraction  =', &
                         optConst(PHASE_PROP_V,RADI_PHASE_CONPART)
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction  =', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_CONPART)
      WRITE(STDOUT,1010) SOLVER_NAME//'   particle diam.=', &
                         optConst(PHASE_PROP_D,RADI_PHASE_CONPART)
    ENDIF
    IF (optConst(PHASE_PROP_V,RADI_PHASE_DISPART) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase     = discrete particles'
      WRITE(STDOUT,1010) SOLVER_NAME//'   vol.fraction  =', &
                         optConst(PHASE_PROP_V,RADI_PHASE_DISPART)
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction .=', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_DISPART)
      WRITE(STDOUT,1010) SOLVER_NAME//'   particle diam.=', &
                         optConst(PHASE_PROP_D,RADI_PHASE_DISPART)
    ENDIF

    extCoef = 0._RFREAL
    DO n = 1,NPHASE
      extCoef = extCoef + optConst(PHASE_PROP_V,n)* &
                optConst(PHASE_PROP_Q,n)/optConst(PHASE_PROP_D,n)
    ENDDO
    extCoef = 1.5_RFREAL*extCoef
    WRITE(STDOUT,1010) SOLVER_NAME//' artif. ext.coef.=', extCoef
    WRITE(STDOUT,1010) SOLVER_NAME//' trshld ext.coef.=', RADI_REAL_ECMIN
    IF (extCoef < RADI_REAL_ECMIN ) &
    WRITE(STDOUT,1030) SOLVER_NAME//' radiation computation has no effect'

  ELSEIF (media == RADI_MEDIA_REAL) THEN
    IF (optConst(PHASE_PROP_V,RADI_PHASE_GAS) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase   = gas'
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction .=', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_GAS)
    ENDIF
    IF (optConst(PHASE_PROP_V,RADI_PHASE_CONPART) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase   = continuum particles'
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction .=', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_CONPART)
    ENDIF
    IF (optConst(PHASE_PROP_V,RADI_PHASE_DISPART) > RADI_REAL_SMALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' media phase   = discrete particles'
      WRITE(STDOUT,1010) SOLVER_NAME//'   Q-extinction .=', &
                         optConst(PHASE_PROP_Q,RADI_PHASE_DISPART)
    ENDIF ! media phase

    WRITE(STDOUT,1010) SOLVER_NAME//' trshld ext.coef.=', RADI_REAL_ECMIN
  ENDIF   ! radiation media

! output selected parameters pertinent to RTE (general) method to monitor

  IF ((radiModel == RADI_MODEL_FLDSRC) .OR. &
      (radiModel == RADI_MODEL_FLDTRAN)) THEN

    IF (fluxLim == FLD_LIM_NONE) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' flux limiter    = none'
    ELSEIF (fluxLim == FLD_LIM_LP) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' flux limiter    = Levermore-Pomraning'
    ENDIF

    IF (radiModel == RADI_MODEL_FLDTRAN) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' FLDTRAN Numeric :'
      WRITE(STDOUT,1020) SOLVER_NAME//' FLD smoocf      =',smoocf
      IF (discrType == FLD_DISCR_CEN) THEN
        IF (discrOrder == FLD_DISCR_ORD1) &
        WRITE(STDOUT,1030) SOLVER_NAME//' FLD spc discr.  = 1st order central'
        IF (discrOrder == FLD_DISCR_ORD2) &
        WRITE(STDOUT,1030) SOLVER_NAME//' FLD spc discr.  = 2nd order central'
        WRITE(STDOUT,1020) SOLVER_NAME//' FLD k2          =',k2
        WRITE(STDOUT,1020) SOLVER_NAME//' FLD 1/k4        =',ook4
      ELSEIF (discrType == FLD_DISCR_UPW) THEN
        IF (discrOrder == FLD_DISCR_ORD1) &
        WRITE(STDOUT,1030) SOLVER_NAME//' FLD spc discr.  = 1st order upwind'
        IF (discrOrder == FLD_DISCR_ORD2) &
        WRITE(STDOUT,1030) SOLVER_NAME//' FLD spc discr.  = 2nd order upwind'
      ENDIF
    ENDIF

  ELSEIF ((radiModel == RADI_MODEL_RTEGRAY) .OR. &
          (radiModel == RADI_MODEL_RTEBAND)) THEN

    IF (solMethod == RADI_NUM_DOM4) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' solution method = discrete ordinate S4'
    ELSEIF (solMethod == RADI_NUM_DOM8) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' solution method = discrete ordinate S8'
    ELSEIF (solMethod == RADI_NUM_DOM16) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' solution method = discrete ordinate S16'
    ELSEIF (solMethod == RADI_NUM_FVM) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' solution method = finite volume method'
    ENDIF

    IF (solMethod == RADI_NUM_FVM) THEN
      WRITE(STDOUT,1020) SOLVER_NAME//' N polar angles  =', nPol
      WRITE(STDOUT,1020) SOLVER_NAME//' N azimuthal     =', nAzi
    ENDIF ! solMethod
  ENDIF   ! radiModel

! general derived parameters

  WRITE(STDOUT,1020)   SOLVER_NAME//' N intens.angles =', nPol
  WRITE(STDOUT,  * )   SOLVER_NAME//'      polar angles    =', &
                                      angles(:,RADI_ANGLE_POLAR)
  WRITE(STDOUT,  * )   SOLVER_NAME//'      polar angles    =', &
                                      angles(:,RADI_ANGLE_AZIMU)

! finish ----------------------------------------------------------------------

  CALL DeregisterFunction( global )

1005 FORMAT(/,4X,A)
1010 FORMAT(6X,A,E12.5)
1020 FORMAT(6X,A,I4)
1030 FORMAT(6X,A)

END SUBROUTINE RADI_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_PrintUserInput.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:11:15  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.7  2004/09/22 01:31:35  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.6  2004/09/18 17:41:27  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.5  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2003/07/30 22:23:41  wasistho
! enter part and smoke data into radiation
!
! Revision 1.3  2003/07/23 03:13:58  wasistho
! cured baby illness
!
! Revision 1.2  2003/07/18 01:39:08  wasistho
! removed bcModel from input data structure
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







