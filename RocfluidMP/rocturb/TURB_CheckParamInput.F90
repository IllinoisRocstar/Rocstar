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
! Purpose: Check TURB parameters either specified by user or set in the code.
!
! Description: The checking includes the existency and order of parameters.
!
! Input: regions = input parameters contained in turbInput of all regions.
!
! Output: Error msg.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_CheckParamInput.F90,v 1.19 2009/08/26 12:28:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CheckParamInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModBndPatch, ONLY   : t_patch
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy
#endif
  USE ModTurbulence
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, m, iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  TYPE(t_mixt_input), POINTER  :: mixtInput  
  TYPE(t_turb_input), POINTER  :: input  
  TYPE(t_patch), POINTER       :: patch1
  LOGICAL :: turbInactive, fixedGrid, moveGrid, wRansActive, noNsWall

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, maxRange
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_CheckParamInput.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_CheckParamInput',&
  'TURB_CheckParamInput.F90' )

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Entering TURB_CheckParamInput...'
  END IF ! global%verbLevel

! check fixed parameters setting ---------------------------------------------

! fixed parameters for general rocturb

  IF ((YCOORD - XCOORD)/=1 .OR. &
      (ZCOORD - YCOORD)/=1) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'XCOORD,YCOORD,ZCOORD' )
  ENDIF

#ifdef RFLO
  IF (NDIR /= 3) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'NDIR /= 3' )
  ENDIF
  IF ((DIRK-DIRJ /= 1).OR.(DIRJ-DIRI /= 1).OR.(DIRI /= 1)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'DIRI,DIRJ,DIRK' )
  ENDIF
#endif

  IF ((CV_TURB_VVEL - CV_TURB_UVEL /= 1) .OR. &
      (CV_TURB_WVEL - CV_TURB_VVEL /= 1)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'CV_TURB_...' )
  ENDIF

  IF ((CV_TURB_XMOM - CV_TURB_DENS /= 1) .OR. &
      (CV_TURB_YMOM - CV_TURB_XMOM /= 1) .OR. &
      (CV_TURB_ZMOM - CV_TURB_YMOM /= 1)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'CV_TURB_...' )
  ENDIF

  IF (((GR_TURB_VX - GR_TURB_UX)/=1).OR.((GR_TURB_WX - GR_TURB_VX)/=1).OR. &
      ((GR_TURB_UY - GR_TURB_WX)/=1).OR.((GR_TURB_VY - GR_TURB_UY)/=1).OR. &
      ((GR_TURB_WY - GR_TURB_VY)/=1).OR.((GR_TURB_UZ - GR_TURB_WY)/=1).OR. &
      ((GR_TURB_VZ - GR_TURB_UZ)/=1).OR.((GR_TURB_WZ - GR_TURB_VZ)/=1)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'GR_TURB_...' )
  ENDIF

#ifdef RFLU
  IF ((CV_TURB_DENS /= CV_MIXT_DENS) .OR. &  ! relevant to LesLij, LesMij, LesHij,
      (CV_TURB_XMOM /= CV_MIXT_XMOM) .OR. &  ! etc
      (CV_TURB_YMOM /= CV_MIXT_YMOM) .OR. &
      (CV_TURB_ZMOM /= CV_MIXT_ZMOM)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'CV_TURB_.../= CV_MIXT_...' )
  ENDIF

!  IF ((CV_TURB_DENS /= 1) .OR. &    ! relevant to RFLU_InterpCells2Faces, etc
!      (E11          /= 1)) THEN
!    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'CV_TURB_DES or E11 /= 1' )
!  ENDIF
#endif

! fixed parameters pertinent to LES

  IF ((FILWIDTH_ONE  /= 1) .OR. &
      (FILWIDTH_TWO  /= 2) .OR. &
      (FILWIDTH_FOUR /= 4)) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'FILWIDTH_...' )
  ENDIF

! fixed parameters pertinent to WLM 

  IF (WLM_VALS_XIX/=1 .OR. WLM_VALS_TAUUX/=10 .OR. &
     (WLM_VALS_ETX-WLM_VALS_XIX)/=1 .OR. (WLM_VALS_ZTX-WLM_VALS_ETX)/=1 .OR. &
     (WLM_VALS_ETY-WLM_VALS_XIY)/=1 .OR. (WLM_VALS_ZTY-WLM_VALS_ETY)/=1 .OR. &
     (WLM_VALS_ETZ-WLM_VALS_XIZ)/=1 .OR. (WLM_VALS_ZTZ-WLM_VALS_ETZ)/=1 .OR. &
     (WLM_VALS_XIY-WLM_VALS_ZTX)/=1 .OR. (WLM_VALS_XIZ-WLM_VALS_ZTY)/=1 .OR. &
     (WLM_VALS_TAUUY-WLM_VALS_TAUUX)/=1.OR.(WLM_VALS_TAUUZ-WLM_VALS_TAUUY)/=1.OR. &
     (WLM_VALS_TAUVY-WLM_VALS_TAUVX)/=1.OR.(WLM_VALS_TAUVZ-WLM_VALS_TAUVY)/=1.OR. &
     (WLM_VALS_TAUWY-WLM_VALS_TAUWX)/=1.OR.(WLM_VALS_TAUWZ-WLM_VALS_TAUWY)/=1.OR. &
     (WLM_VALS_TAUVX-WLM_VALS_TAUUZ)/=1.OR.(WLM_VALS_TAUWX-WLM_VALS_TAUVZ)/=1) THEN
    CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__,'WLM_VALS_...' )
  ENDIF

! check turbulence parameters selection

  turbInactive = .FALSE.
  fixedGrid    = .FALSE.
  moveGrid     = .FALSE.
  wRansActive  = .FALSE.
  noNsWall     = .TRUE.

#ifdef RFLO
  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif

      mixtInput => regions(iReg)%mixtInput
      input     => regions(iReg)%turbInput

      IF ((mixtInput%turbModel /= TURB_MODEL_NONE).AND. &
          (mixtInput%turbModel /= TURB_MODEL_FIXSMAG).AND. &
          (mixtInput%turbModel /= TURB_MODEL_SCALSIM).AND. &
          (mixtInput%turbModel /= TURB_MODEL_DYNSMAG).AND. &
          (mixtInput%turbModel /= TURB_MODEL_DYNMIXD).AND. &
          (mixtInput%turbModel /= TURB_MODEL_SA)     .AND. &
          (mixtInput%turbModel /= TURB_MODEL_DESSA)  .AND. &
          (mixtInput%turbModel /= TURB_MODEL_HDESSA)) THEN
         CALL ErrorStop( global,ERR_TURB_MODEL,__LINE__ )
      ENDIF
#ifdef RFLU
      IF (mixtInput%turbModel==TURB_MODEL_DYNMIXD) THEN
        CALL ErrorStop( global,ERR_TURB_MODEL,__LINE__, &
          'LES Dynamic Mixed model is not available with unstructured grid' )
      ENDIF
#endif
      IF ((input%modelClass /= MODEL_RANS) .AND. (input%nCv > 0)) THEN
        CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__, &
             'number of conservative variables > 0 only for RANS/DES' )
      ENDIF

      IF ((input%calcVort /= CALCVORT_NO) .AND. &
          (input%calcVort /= CALCVORT_FDT) .AND. &
          (input%calcVort /= CALCVORT_SDT)) THEN 
        CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__,'CALCVORTIC: 0, 1 or 2' )
      ENDIF

#ifdef GENX
      IF (input%calcVort == CALCVORT_NO) THEN
        CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
        'Vorticities should always be computed in Genx, CALCVORTIC cannot be 0' )
      ENDIF
#endif

      IF ((input%nOutField > 1).AND.(input%calcVort == CALCVORT_NO)) THEN
        CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
                        'CALCVORTIC should be > 0 for OUTPUTNUMBER > 1' )
      ENDIF

      IF (input%nZof > ZOF_NELM) THEN
        CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__, &
                        'ZOF_NELM < input%nZof, increase the former' )
      ENDIF

! --- pertinent to LES

      IF (input%modelClass == MODEL_LES) THEN
        IF ((mixtInput%turbModel==TURB_MODEL_FIXSMAG).OR. &
            (mixtInput%turbModel==TURB_MODEL_SCALSIM).OR. &
            (mixtInput%turbModel==TURB_MODEL_DYNSMAG).OR. &
            (mixtInput%turbModel==TURB_MODEL_DYNMIXD)) THEN
        ELSE
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
            'selected turbulence model is not of LES class' )
        ENDIF

        IF ((input%nOutField < 1).OR.(input%nOutField > MAXOUTFLD_LES)) THEN
          CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
                          'OUTPUTNUMBER for LES out of range' )
        ENDIF

        IF ((input%cSmag < 0._RFREAL).OR.(input%cSmag > 0.2_RFREAL)) THEN
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                          '0<= CSMAGORINSKY <=0.2' )
        ENDIF

#ifdef RFLO
        IF ((input%filterType /= FILTYPE_UNIFORM) .AND. &
            (input%filterType /= FILTYPE_NONUNIF)) THEN
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                          'FILTERTYPE: 0 or 1' )
        ENDIF

        IF ((input%deltaType /= DELTYPE_CBRT) .AND. &
            (input%deltaType /= DELTYPE_SQRT)) THEN
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                          'DELTATYPE: 0 or 1' )
        ENDIF

        DO m = DIRI,DIRK
          IF ((input%filterWidth(m) /= FILWIDTH_ZERO) .AND. &
              (input%filterWidth(m) /= FILWIDTH_ONE)  .AND. &
              (input%filterWidth(m) /= FILWIDTH_TWO)) THEN
            CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                            'FILTERWIDTH: 0,1 or 2' )
          ENDIF
          IF ((input%homDir(m) /= OFF) .AND. &
              (input%homDir(m) /= ACTIVE)) THEN
            CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                            'HOMOGENDIR: 0 or 1' )
          ENDIF
        ENDDO
#endif
#ifdef RFLU
        IF ((input%filterWidth(DIRI) /= FILWIDTH_ZERO) .AND. &
            (input%filterWidth(DIRI) /= FILWIDTH_ONE)  .AND. &
            (input%filterWidth(DIRI) /= FILWIDTH_TWO)) THEN
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                          'FILTERWIDTH: 0,1 or 2' )
        ENDIF
#endif
        IF ((input%engModel /= OFF) .AND. &
            (input%engModel /= ACTIVE)) THEN
          CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__, &
                          'ENERGYMODEL: 0 or 1' )
        ENDIF

      ENDIF  ! modelClass LES

! --- pertinent to RaNS/DES

      IF (input%modelClass == MODEL_RANS) THEN
        IF ((mixtInput%turbModel==TURB_MODEL_SA).OR. &
            (mixtInput%turbModel==TURB_MODEL_DESSA).OR. &
            (mixtInput%turbModel==TURB_MODEL_HDESSA)) THEN
          IF (input%wDistMethod /= WDIST_DIRECT) THEN
            CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
              'only direct WALLDISTMETHOD (0) currently valid for SA/DES' )
          ENDIF
!          IF ((mixtInput%moveGrid .eqv. .true.) .AND. (input%wDistMethod==WDIST_DIRECT)) THEN
!            CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
!              'direct WALLDISTMETHOD (0) is not efficient for moving grid' )
!          ENDIF
          IF (input%functV1 /= SA_FV1_POW3 .AND. &
              input%functV1 /= SA_FV1_POW2) THEN
            CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
              'selected formula for SA function fv1 is invalid' )
          ENDIF
        ELSE
          CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
            'selected turbulence model is not of RANS/DES class' )
        ENDIF

        IF ((input%nOutField < 1).OR.(input%nOutField > MAXOUTFLD_RANS)) THEN
          CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
                          'OUTPUTNUMBER for RANS out of range' )
        ENDIF

#ifdef RFLO
        IF (input%spaceDiscr /= RANS_DISCR_CEN .AND. &
            input%spaceDiscr /= RANS_DISCR_UPW) THEN
          CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
            'selected RaNS space discretization is not defined' )
        ENDIF
#endif
#ifdef RFLU
        IF (input%spaceDiscr /= RANS_DISCR_UPW) THEN
          CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
            'selected RaNS space discretization is not defined in Rocflu' )
        ENDIF
#endif

        IF (input%spaceOrder /= RANS_DISCR_ORD1 .AND. &
            input%spaceOrder /= RANS_DISCR_ORD2) THEN
          CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
            'selected RaNS space discretization order is not defined' )
        ENDIF
      ENDIF

! --- pertinent to WLM
 
#ifdef RFLO
      iLev = 1  ! check input based on finest level
      CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                jdcbeg,jdcend,kdcbeg,kdcend )

      DO iPatch=1,regions(iReg)%nPatches
        patch1 => regions(iReg)%levels(iLev)%patches(iPatch)
#endif
#ifdef RFLU
      DO iPatch=1,regions(iReg)%grid%nPatches
        patch1 => regions(iReg)%patches(iPatch)
#endif
        IF (patch1%bcType>=BC_NOSLIPWALL .AND. &
            patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
          
          IF (patch1%valBola%nSwitches <= 0) THEN
            CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__, &
                           'nSwitches in bcvalues type valBola should be > 0' )
          ENDIF
          
          IF (patch1%valBola%switches(WLM_INPUT_MODEL) < 0 .OR. &
              patch1%valBola%switches(WLM_INPUT_MODEL) > WLM_MODEL_EXTERN) THEN
            CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__,'MODEL: 0-3' )
          ENDIF
          
          IF (patch1%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_EXTERN) THEN
            CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
            'Ext. tau_wall not ready yet' )
          ENDIF

          IF (patch1%valBola%switches(WLM_INPUT_MODEL)/=WLM_MODEL_NOMODEL) THEN

            IF (MINVAL(patch1%valBola%vals(:,WLM_VALS_ROUGH)) < 0._RFREAL) THEN
              CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                           'ROUGHNESS: > 0.0' )
            ENDIF
#ifdef RFLO
            IF (patch1%valBola%switches(WLM_INPUT_REFPOINT) < 1) THEN
              CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                             'WLM reference point should be >= 1' )
            ENDIF
            maxRange = idcend-idcbeg+ 1
            IF ((patch1%lbound == 1 .OR. patch1%lbound == 2) .AND. &
                 patch1%valBola%switches(WLM_INPUT_REFPOINT) > maxRange) THEN
              CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                             'WLM reference point exceeds region max. range' )
            ENDIF
            maxRange = jdcend-jdcbeg+ 1
            IF ((patch1%lbound == 3 .OR. patch1%lbound == 4) .AND. &
                 patch1%valBola%switches(WLM_INPUT_REFPOINT) > maxRange) THEN
              CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                             'WLM reference point exceeds region max. range' )
            ENDIF
            maxRange = kdcend-kdcbeg+ 1
            IF ((patch1%lbound == 5 .OR. patch1%lbound == 6) .AND. &
                 patch1%valBola%switches(WLM_INPUT_REFPOINT) > maxRange) THEN
              CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                             'WLM reference point exceeds region max. range' )
            ENDIF
#endif
#ifdef RFLU
            CALL ErrorStop( global,ERR_TURB_WLMINPUT,__LINE__, &
                           'Wall Layer Model is not available yet in Rocflu' )
#endif
          ENDIF   ! switch
        ENDIF     ! bcType
      ENDDO       ! iPatch    

! --- check turbulence statistics input

      IF (input%nSt > ST_TURB_NVAR) THEN
        CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__, &
                       'increase ST_TURB_NVAR in parameters module file' )
      ENDIF

#ifdef STATS
      IF (global%turbNStat > (input%nFixSt+input%nSv+input%nSt)) THEN
        CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
            'TURBNSTAT larger than allocated; nSv (stress), nSt (stats) may be 0' )
      ENDIF

      IF ((global%doStat==1) .AND. &
          (mixtInput%turbModel/=TURB_MODEL_DYNSMAG) .AND. &
          (mixtInput%turbModel/=TURB_MODEL_DYNMIXD)) THEN
        DO m = 1,global%turbNStat
          IF (global%turbStatId(1,m)==0 .AND. global%turbStatId(2,m)>2 ) THEN
            CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
            'TURBSTATID > 02 only for dyn.LES; chg setting if otherwise desired' )
          ENDIF
        ENDDO
      ENDIF
#endif

#ifdef RFLO
    ENDIF ! region active
#endif
  ENDDO   ! iReg

! global checking --------------------------------------------------------

#ifdef RFLO
  DO iReg=1,global%nRegions    ! serial process
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif
    mixtInput => regions(iReg)%mixtInput
    input     => regions(iReg)%turbInput

! - general: turbInactive and moveGrid vs fixedGrid

    IF (mixtInput%turbModel == TURB_MODEL_NONE) THEN
      turbInactive = .TRUE.
    ENDIF
    IF       (mixtInput%moveGrid .eqv. .true.)  moveGrid = .TRUE.
    IF (.NOT. (mixtInput%moveGrid .eqv. .true.)) fixedGrid = .TRUE.

! - RaNS/DES: wRansActive and noNsWall

    IF (input%modelClass == MODEL_RANS) THEN
      IF ((mixtInput%turbModel==TURB_MODEL_SA).OR. &
          (mixtInput%turbModel==TURB_MODEL_DESSA).OR. &
          (mixtInput%turbModel==TURB_MODEL_HDESSA)) THEN
        wRansActive = .TRUE.          
      ENDIF
    ENDIF

# ifdef RFLO
    DO iPatch=1,regions(iReg)%nPatches
      patch1 => regions(iReg)%levels(iLev)%patches(iPatch)
#endif
# ifdef RFLU
    DO iPatch=1,regions(iReg)%grid%nPatches
      patch1 => regions(iReg)%patches(iPatch)
#endif
      IF (patch1%bcType>=BC_NOSLIPWALL .AND. &
          patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
        noNsWall = .FALSE.
      ENDIF
    ENDDO

  ENDDO ! iReg

#ifdef GENX
  IF (turbInactive .eqv. .true.) THEN
    CALL ErrorStop( global,ERR_TURB_REGION,__LINE__, &
      'For Genx, turbulence should be active in all regions or none at all' )
  ENDIF
#endif
  IF ((fixedGrid .eqv. .true.) .AND. (moveGrid .eqv. .true.)) THEN
    CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
      'QUESTION: Is it safe to only update metrics in regions that move?' )
  ENDIF

  IF ((wRansActive .eqv. .true.) .AND. (noNsWall .eqv. .true.)) THEN
    CALL ErrorStop( global,ERR_TURB_INPUT,__LINE__, &
      'RaNS model selected needs wall distance but no ns-wall presents' )
  ENDIF

! finalize -------------------------------------------------------------------

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Leaving TURB_CheckParamInput.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CheckParamInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_CheckParamInput.F90,v $
! Revision 1.19  2009/08/26 12:28:52  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.18  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2006/02/04 05:00:18  wasistho
! added enter and leave statements
!
! Revision 1.15  2006/01/12 09:50:07  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.14  2005/12/29 19:47:14  wasistho
! bug fixed start/end index of ireg for rflu
!
! Revision 1.13  2005/03/09 06:34:45  wasistho
! incorporated HDESSA
!
! Revision 1.12  2005/03/08 01:54:56  wasistho
! fixed bug, added if dostat cheking in turb-statistics
!
! Revision 1.11  2004/12/10 01:03:13  wasistho
! trap turbStatId error in connection with turbmodel selected
!
! Revision 1.10  2004/08/12 20:58:39  wasistho
! check roughness value only when WLM is active
!
! Revision 1.9  2004/08/07 01:08:03  wasistho
! error msg if Dyn.Mixed selected in Rocflu
!
! Revision 1.8  2004/07/03 01:44:00  wasistho
! removed restriction to only fixed Smag in Rocflu
!
! Revision 1.7  2004/05/28 01:59:43  wasistho
! update unstructured grid LES
!
! Revision 1.6  2004/04/08 04:00:27  wasistho
! allow DES with moving grid
!
! Revision 1.5  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.4  2004/03/24 03:37:02  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/19 02:45:12  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 21:08:02  wasistho
! changed nomenclature
!
! Revision 1.18  2004/02/26 21:25:50  wasistho
! delete energy sgs warning
!
! Revision 1.17  2004/02/24 21:03:23  wasistho
! modified the warning previously checked in
!
! Revision 1.16  2004/02/24 02:52:18  wasistho
! added warning for unbalance load MPI if LES energy model is used
!
! Revision 1.15  2004/02/19 04:03:42  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.14  2004/02/14 03:43:07  wasistho
! added new WLM parameter: reference point
!
! Revision 1.13  2004/02/11 03:24:40  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.12  2003/12/24 00:04:36  wasistho
! modify WLM checking to satisfy Turing
!
! Revision 1.11  2003/10/26 00:56:14  wasistho
! bug fixed
!
! Revision 1.10  2003/10/26 00:09:26  wasistho
! added multiple discr.types and order
!
! Revision 1.9  2003/10/07 02:05:04  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.8  2003/09/11 03:07:30  wasistho
! removed ASSOCIATED(patch1%valBola%switches)
!
! Revision 1.7  2003/08/29 20:30:35  wasistho
! Added check for vorticity computation in Genx
!
! Revision 1.6  2003/08/06 15:56:05  wasistho
! added vorticities computation
!
! Revision 1.5  2003/07/18 20:12:56  wasistho
! put ifdef STATS
!
! Revision 1.4  2003/07/17 01:23:10  wasistho
! prepared rocturb for Genx
!
! Revision 1.3  2003/05/31 01:46:21  wasistho
! installed turb. wall layer model
!
! Revision 1.2  2003/05/24 02:11:05  wasistho
! turbulence statistics expanded
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







