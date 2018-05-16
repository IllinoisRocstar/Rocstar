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
! Purpose: Check RADI parameters either specified by user or set in the code.
!
! Description: Fixed parameters and relevancy of physical parameters are
!              checked.
!
! Input: regions = input parameters contained in radiInput of all regions.
!
! Output: Error msg for inconsistency.
!
! Notes: If radiation is used in Genx, it should be active in all regions
!        or not at all. When this routine is called, radiation is active
!        at least in one region.
!
!******************************************************************************
!
! $Id: RADI_CheckParamInput.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

#ifdef RFLO
SUBROUTINE RADI_CheckParamInput( regions )
#endif
#ifdef RFLU
SUBROUTINE RADI_CheckParamInput
#endif

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModRadiation, ONLY  : t_radi_input
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER      :: global
  TYPE(t_radi_input), POINTER  :: input

  LOGICAL :: radiUnused

  REAL(RFREAL), POINTER :: optConst(:,:)

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RADI_CheckParamInput',&
  'RADI_CheckParamInput.F90' )

! check fixed parameters setting ---------------------------------------------

  IF ((YCOORD - XCOORD)/=1 .OR. &
      (ZCOORD - YCOORD)/=1) THEN
    CALL ErrorStop( global,ERR_RADI_FIXPARAM,__LINE__,'XCOORD,YCOORD,ZCOORD' )
  ENDIF

#ifdef RFLO
  IF ((KCOORD-JCOORD /= 1).OR.(JCOORD-ICOORD /= 1).OR.(ICOORD /= 1)) THEN
    CALL ErrorStop( global,ERR_RADI_FIXPARAM,__LINE__,'I,J,KCOORD' )
  ENDIF
#endif

  IF (RADI_COEFF_EXTINCT /=1 .OR. &
      RADI_COEFF_SCATTER /=2 .OR. &
      RADI_COEFF_NCOMP /=2) THEN
    CALL ErrorStop( global,ERR_RADI_FIXPARAM,__LINE__,'RADI_COEFF_...' )
  ENDIF

  IF (RADI_ANGLE_NCOMP /=2) THEN
    CALL ErrorStop( global,ERR_RADI_FIXPARAM,__LINE__,'RADI_ANGLE_...' )
  ENDIF

! check RADI parameter selection regionwise ---------------------------------

  radiUnused = .false.

#ifdef RFLO
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

! --- radi input check

      input    => regions(iReg)%radiInput
      optConst => input%optConst

! --- radiation model

      IF ((input%radiModel /= RADI_MODEL_NONE)    .AND. &
          (input%radiModel /= RADI_MODEL_ROSS)    .AND. &
          (input%radiModel /= RADI_MODEL_FLDSRC)  .AND. &
          (input%radiModel /= RADI_MODEL_FLDTRAN) .AND. &
          (input%radiModel /= RADI_MODEL_RTEGRAY) .AND. &
          (input%radiModel /= RADI_MODEL_RTEBAND)) THEN
        CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
             'radiation model selected does not exist' )
      ENDIF

      IF (regions(iReg)%mixtInput%flowModel /= FLOW_NAVST) THEN
        IF ((input%radiModel == RADI_MODEL_ROSS)   .OR. &
            (input%radiModel == RADI_MODEL_FLDSRC) .OR. &
            (input%radiModel == RADI_MODEL_FLDTRAN)) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
              'diffusion approximation models currently only run with NS' )
        ENDIF
      ENDIF

! --- radiation active

      IF (input%radiModel /= RADI_MODEL_NONE) THEN

! ----- active radiation model
        IF ((input%radiModel == RADI_MODEL_FLDTRAN) .OR. &
            (input%radiModel == RADI_MODEL_RTEGRAY) .OR. &
            (input%radiModel == RADI_MODEL_RTEBAND)) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
               'radiation model selected is not ready yet' )
        ENDIF

! ----- check for real media, consistency with rocpart and rocsmoke
        IF (input%media == RADI_MEDIA_REAL) THEN
          IF (optConst(PHASE_PROP_V,RADI_PHASE_DISPART) >RADI_REAL_SMALL) THEN
            IF (.NOT. global%plagUsed) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
              'real media contains discrete particles but rocpart is off' )
            ENDIF
          ENDIF   ! optConst
          IF (optConst(PHASE_PROP_V,RADI_PHASE_CONPART) >RADI_REAL_SMALL) THEN
            IF (.NOT. global%peulUsed) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
              'real media contains continuum particles but rocsmoke is off' )
            ENDIF
          ENDIF
        ENDIF     ! media

! ----- check consistency of optical constants, independent of media
        IF (optConst(PHASE_PROP_V,RADI_PHASE_GAS) <RADI_REAL_SMALL .OR. &
            optConst(PHASE_PROP_Q,RADI_PHASE_GAS) <RADI_REAL_SMALL) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
               'Gas: ext.coef. and vol.frac. should be > 0 or let be default' )
        ENDIF
        IF (optConst(PHASE_PROP_V,RADI_PHASE_CONPART) >RADI_REAL_SMALL .AND. &
            optConst(PHASE_PROP_Q,RADI_PHASE_CONPART) <RADI_REAL_SMALL) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
               'Con. particles: extinction efficiency should be higher' )
        ENDIF
        IF (optConst(PHASE_PROP_V,RADI_PHASE_DISPART) >RADI_REAL_SMALL .AND. &
            optConst(PHASE_PROP_Q,RADI_PHASE_DISPART) <RADI_REAL_SMALL) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
               'Disc. particles: extinction efficiency should be higher' )
        ENDIF
        IF (optConst(PHASE_PROP_D,RADI_PHASE_GAS    ) <RADI_REAL_SMALL .OR. &
            optConst(PHASE_PROP_D,RADI_PHASE_DISPART) <RADI_REAL_SMALL .OR. &
            optConst(PHASE_PROP_D,RADI_PHASE_CONPART) <RADI_REAL_SMALL) THEN
        CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
             'gas molecule and particle diameters should be > (machine) zero' )
        ENDIF

! ----- intensity angles
        IF ((LBOUND(input%angles,1) /= 1) .OR. &
            (UBOUND(input%angles,1) /= input%nAng)) THEN
          CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
               'number of intensity angles is inconsistent' )
        ENDIF

! ----- radiation numerical solution method of RTE radiation model
        IF (input%radiModel == RADI_MODEL_FLDTRAN) THEN
#ifdef RFLO
          IF (input%spaceDiscr /= FLD_DISCR_CEN) THEN
            CALL ErrorStop( global,ERR_RADI_FLDINPUT,__LINE__, &
              'selected FLD space discretization is not defined' )
          ENDIF
#endif
          IF (input%spaceOrder /= FLD_DISCR_ORD2) THEN
            CALL ErrorStop( global,ERR_RADI_FLDINPUT,__LINE__, &
              'selected FLD space discretization order is not defined' )
          ENDIF
        ENDIF

! ----- radiation numerical solution method of RTE radiation model
        IF ((input%radiModel == RADI_MODEL_RTEGRAY) .OR. &
            (input%radiModel == RADI_MODEL_RTEBAND)) THEN

! ------- DOM:
          IF ((input%solMethod /= RADI_NUM_DOM4)  .AND. &
              (input%solMethod /= RADI_NUM_DOM8)  .AND. &
              (input%solMethod /= RADI_NUM_DOM16) .AND. &
              (input%solMethod /= RADI_NUM_FVM))  THEN
            CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                 'radiation solution method selected does not exist' )
          ENDIF
          IF ((input%solMethod == RADI_NUM_DOM4)  .OR. &
              (input%solMethod == RADI_NUM_DOM8)  .OR. &
              (input%solMethod == RADI_NUM_DOM16) .OR. &
              (input%solMethod == RADI_NUM_FVM))  THEN
            CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                 'radiation solution method selected is not ready yet' )
          ENDIF
          IF ((input%solMethod == RADI_NUM_DOM4) .AND. &
              (input%nOrdin /= 4)) THEN
            CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                 'number of ordinates inconsistence with method selected' )
          ENDIF
          IF ((input%solMethod == RADI_NUM_DOM8) .AND. &
              (input%nOrdin /= 8)) THEN
            CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                 'number of ordinates inconsistence with method selected' )
          ENDIF
          IF ((input%solMethod == RADI_NUM_DOM16) .AND. &
              (input%nOrdin /= 16)) THEN
            CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                 'number of ordinates inconsistence with method selected' )
          ENDIF
          IF ((input%solMethod == RADI_NUM_DOM4)  .OR. &
              (input%solMethod == RADI_NUM_DOM8)  .OR. &
              (input%solMethod == RADI_NUM_DOM16)) THEN
            IF (input%nAng /= input%nOrdin) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                   'DOM: nAngle /= nOrdinate' )
            ENDIF
          ENDIF

! ------- FVM:
          IF (input%solMethod == RADI_NUM_FVM) THEN
            IF (input%nPol < 2 .OR. input%nPol > 50) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                   'number of polar angles is out of range [2-50]' )
            ENDIF
            IF (input%nAzi < 4 .OR. input%nAzi > 100) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                   'number of azimuthal angles is out of range [4-100]' )
            ENDIF
            IF (input%nAng /= (input%nPol+1)*(input%nAzi+1)) THEN
              CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
                   'FVM: nAngle /= (nPolar+1)*(nAzimuthal+1)' )
            ENDIF
          ENDIF  ! solMethod
        ENDIF    ! radiModel RTE

      ENDIF      ! radiModel active

! --- assign value to radiUnused
      IF (.NOT. regions(iReg)%mixtInput%radiUsed) THEN
        radiUnused = .true.
      ENDIF
#endif
#ifdef RFLU
      input => radiInput
      IF (input%radiModel /= RADI_MODEL_NONE) THEN
        CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
             'RFLU-RADI is not ready yet' )
      ENDIF
#endif

#ifdef RFLO
    ENDIF ! region active
  ENDDO   ! iReg
#endif

#ifdef GENX
  IF (radiUnused) THEN
    CALL ErrorStop( global,ERR_RADI_INPUT,__LINE__, &
         'For Genx, radiation should be active in all regions or none at all' )
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_CheckParamInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_CheckParamInput.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:10:50  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.8  2004/09/22 01:30:33  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.7  2004/09/18 17:40:44  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.6  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2003/08/07 21:58:23  wasistho
! diffusion approx. currently only run with NS
!
! Revision 1.4  2003/07/30 22:22:31  wasistho
! enter part and smoke data into radiation
!
! Revision 1.3  2003/07/23 03:13:43  wasistho
! cured baby illness
!
! Revision 1.2  2003/07/18 01:44:18  wasistho
! check consistency with rocpart and rocsmoke
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







