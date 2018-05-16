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
! Purpose: Write out TURB user input for checking purposes.
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
! $Id: TURB_PrintUserInput.F90,v 1.9 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_PrintUserInput( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER       :: turbModel, modelClass
  INTEGER       :: calcVort

! ... WLM variables
  INTEGER       :: wallModel, wlmRefPoint
  REAL(RFREAL)  :: wallRough

! ... LES variables
  INTEGER       :: engModel, filterWidth(DIRI:DIRK)
  REAL(RFREAL)  :: cSmag, xyzSmag(XCOORD:ZCOORD)
#ifdef RFLO
  INTEGER       :: filterType, deltaType, homDir(DIRI:DIRK)
#endif

! ... RaNS/DES variables
  INTEGER       :: wDistMeth, functV1, discrType, discrOrder
  REAL(RFREAL)  :: cb1, cb2, cw2, cw3, cv1, rsigma, rkappa, cw1, cdes
  REAL(RFREAL)  :: smoocf, k2, ook4

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_PrintUserInput.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_PrintUserInput',&
  'TURB_PrintUserInput.F90' )

! get pointers ---------------------------------------------------------------

  turbModel      = region%mixtInput%turbModel
  calcVort       = region%turbInput%calcVort

! WLM
  wallModel      = region%turbInput%wallModel
  wlmRefPoint    = region%turbInput%wlmRefPoint
  wallRough      = region%turbInput%wallRough

! LES
  cSmag          = region%turbInput%cSmag
  xyzSmag(:)     = region%turbInput%xyzSmag(:)
  modelClass     = region%turbInput%modelClass
  filterWidth(:) = region%turbInput%filterWidth(:)
  engModel       = region%turbInput%engModel
#ifdef RFLO
  filterType     = region%turbInput%filterType
  deltaType      = region%turbInput%deltaType
  homDir(:)      = region%turbInput%homDir(:)
#endif

! RaNS/DES
  IF ((turbModel==TURB_MODEL_SA).OR.(turbModel==TURB_MODEL_DESSA).OR. &
      (turbModel==TURB_MODEL_HDESSA)) THEN
    cb1          = region%turbInput%const(MC_SA_CB1)
    cb2          = region%turbInput%const(MC_SA_CB2)
    cw2          = region%turbInput%const(MC_SA_CW2)
    cw3          = region%turbInput%const(MC_SA_CW3)
    cv1          = region%turbInput%const(MC_SA_CV1)
    rsigma       = region%turbInput%const(MC_SA_RSIG)
    rkappa       = region%turbInput%const(MC_SA_RKAP)
    cw1          = region%turbInput%const(MC_SA_CW1)
    cDes         = region%turbInput%cdes
  ENDIF
  wDistMeth      = region%turbInput%wDistMethod
  functV1        = region%turbInput%functV1
  smoocf         = region%turbInput%smoocf
  discrType      = region%turbInput%spaceDiscr
  k2             = region%turbInput%vis2
  ook4           = 1._RFREAL/region%turbInput%vis4
  discrOrder     = region%turbInput%spaceOrder

! turbulence model selection

  IF (turbModel==TURB_MODEL_FIXSMAG) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = basic Smagorinsky'
  ELSEIF (turbModel==TURB_MODEL_SCALSIM) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = scale similarity'
  ELSEIF (turbModel==TURB_MODEL_DYNSMAG) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = dynamic Smagorinsky'
  ELSEIF (turbModel==TURB_MODEL_DYNMIXD) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = dynamic mixed'
  ELSEIF (turbModel==TURB_MODEL_SA) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = Spalart-Allmaras'
  ELSEIF (turbModel==TURB_MODEL_DESSA) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = DES-SA'
  ELSEIF (turbModel==TURB_MODEL_HDESSA) THEN
    WRITE(STDOUT,1005) SOLVER_NAME//'     turbulence model = Hybrid DES-SA'
  ENDIF
  WRITE(STDOUT,1030) SOLVER_NAME//'     relevant input-data based on above model:'

! model parameters selection

  IF (modelClass==MODEL_LES) THEN

! - LES part 

    WRITE(STDOUT,1030) SOLVER_NAME//'     class of model   = LES'
    IF (engModel==ACTIVE) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'     energy SGS model = active'
    ELSE
      WRITE(STDOUT,1030) SOLVER_NAME//'     energy SGS model = off'
    ENDIF

    IF (turbModel==TURB_MODEL_FIXSMAG) THEN
      WRITE(STDOUT,1020) SOLVER_NAME//'     C-Smagorinsky    =',cSmag
    ENDIF
    IF (turbModel==TURB_MODEL_FIXSMAG .OR. &
        turbModel==TURB_MODEL_DYNSMAG .OR. &
        turbModel==TURB_MODEL_DYNMIXD) THEN
      IF (xyzSmag(XCOORD) > -HUGE(1._RFREAL )/1000._RFREAL) THEN
        WRITE(STDOUT,1020) SOLVER_NAME//'     x-tripping       =',xyzSmag(XCOORD)
      ENDIF
      IF (xyzSmag(YCOORD) > -HUGE(1._RFREAL )/1000._RFREAL) THEN
        WRITE(STDOUT,1020) SOLVER_NAME//'     y-tripping       =',xyzSmag(YCOORD)
      ENDIF
      IF (xyzSmag(ZCOORD) > -HUGE(1._RFREAL )/1000._RFREAL) THEN
        WRITE(STDOUT,1020) SOLVER_NAME//'     z-tripping       =',xyzSmag(ZCOORD)
      ENDIF
    ENDIF
    IF (turbModel/=TURB_MODEL_FIXSMAG) THEN
#ifdef RFLO
      IF (filterType==FILTYPE_UNIFORM) THEN
        WRITE(STDOUT,1030) SOLVER_NAME//'     filter type      = uniform'
      ELSE
        WRITE(STDOUT,1030) SOLVER_NAME//'     filter type      = non-uniform'
      ENDIF
      IF (deltaType==DELTYPE_CBRT) THEN
        WRITE(STDOUT,1030) SOLVER_NAME//'     delta type       = cuberoot formula'
      ELSE
        WRITE(STDOUT,1030) SOLVER_NAME//'     delta type       = sqrroot formula'
      ENDIF

      WRITE(STDOUT,1010) SOLVER_NAME//'     ijk-filter width =', &
      filterWidth(DIRI), filterWidth(DIRJ), &
      filterWidth(DIRK),' grid-spacing'

      WRITE(STDOUT,1010) SOLVER_NAME//'     ijk-homogen dir. =', &
      homDir(DIRI),homDir(DIRJ),homDir(DIRK), &
      ' (1 = homogeneous, 0 = non-homogeneous)'
#endif
#ifdef RFLU
      WRITE(STDOUT,1008) SOLVER_NAME//'     filter width     =', &
      filterWidth(DIRI),' grid-spacing'
#endif
    ENDIF ! turbModel
  ELSE

! - RANS/DES part 

    WRITE(STDOUT,1030) SOLVER_NAME//'     class of model   = RANS or DES'
    IF ((turbModel==TURB_MODEL_SA).OR.(turbModel==TURB_MODEL_DESSA).OR. &
        (turbModel==TURB_MODEL_HDESSA)) THEN
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cb1 =',cb1
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cb2 =',cb2
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cw2 =',cw2
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cw3 =',cw3
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cv1 =',cv1
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. SIG =',1._RFREAL/rsigma
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. KAP =',1._RFREAL/rkappa
      WRITE(STDOUT,1020) SOLVER_NAME//'     model const. cw1 =',cw1
      IF (turbModel==TURB_MODEL_DESSA .OR. turbModel==TURB_MODEL_HDESSA) &
      WRITE(STDOUT,1020) SOLVER_NAME//'     DES spacing coef =',cDes
      IF (wDistMeth == WDIST_DIRECT) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     wall dist.method = direct'
      IF (wDistMeth == WDIST_HIERAR) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     wall dist.method = hierarchical'
      IF (functV1 == SA_FV1_POW3) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     visc. funct. fv1 = power 3 formula'
      IF (functV1 == SA_FV1_POW2) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     visc. funct. fv1 = power 2 formula'
    ENDIF
    WRITE(STDOUT,1030) SOLVER_NAME//'     RANS/DES Numeric :'
    WRITE(STDOUT,1020) SOLVER_NAME//'     RANS smoocf      =',smoocf
    IF (discrType == RANS_DISCR_CEN) THEN
      IF (discrOrder == RANS_DISCR_ORD1) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     RANS spc discr.  = 1st order central'
      IF (discrOrder == RANS_DISCR_ORD2) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     RANS spc discr.  = 2nd order central'
      WRITE(STDOUT,1020) SOLVER_NAME//'     RANS k2          =',k2
      WRITE(STDOUT,1020) SOLVER_NAME//'     RANS 1/k4        =',ook4
    ELSEIF (discrType == RANS_DISCR_UPW) THEN
      IF (discrOrder == RANS_DISCR_ORD1) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     RANS spc discr.  = 1st order upwind'
      IF (discrOrder == RANS_DISCR_ORD2) &
      WRITE(STDOUT,1030) SOLVER_NAME//'     RANS spc discr.  = 2nd order upwind'
    ENDIF
  ENDIF

! vorticity components

  IF (calcVort == CALCVORT_FDT) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'     vorticity sample : per fluid dt'
  ELSEIF (calcVort == CALCVORT_SDT) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'     vorticity sample : per system dt (Genx)'
  ENDIF

! wall model selection

  IF (wallModel==WLM_MODEL_LOGLAY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'     best wall model  = log layer model'
  ELSEIF (wallModel==WLM_MODEL_BNDLAY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'     best wall model  = boundary layer model'
  ELSEIF (wallModel==WLM_MODEL_EXTERN) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'     best wall model  = prescribed wall stress'
  ENDIF
  IF (wallModel/=WLM_MODEL_NOMODEL) THEN
    WRITE(STDOUT,1020) SOLVER_NAME//'     max. ref. point  =', wlmRefPoint
    IF (wallRough > REAL_SMALL) THEN
      WRITE(STDOUT,1020) SOLVER_NAME//'     max. roughness   =', wallRough
    ENDIF
  ENDIF

! finish ----------------------------------------------------------------------

  CALL DeregisterFunction( global )

1005 FORMAT(/,A)
1008 FORMAT(A,I8,A)
1010 FORMAT(A,3I4,A)
1015 FORMAT(A,3I4)
1020 FORMAT(A,E12.5)
1030 FORMAT(A)

END SUBROUTINE TURB_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_PrintUserInput.F90,v $
! Revision 1.9  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/01/18 06:38:46  wasistho
! screen printed tripping locations for eddyvis LES
!
! Revision 1.6  2006/01/12 09:50:15  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.5  2005/03/09 06:35:10  wasistho
! incorporated HDESSA
!
! Revision 1.4  2004/04/08 03:18:03  wasistho
! left justification of screen output
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/19 02:47:33  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.12  2004/02/19 04:03:50  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.11  2004/02/14 03:43:15  wasistho
! added new WLM parameter: reference point
!
! Revision 1.10  2003/10/27 23:13:19  wasistho
! bug fixed in the output of cDes
!
! Revision 1.9  2003/10/26 00:56:24  wasistho
! bug fixed
!
! Revision 1.8  2003/10/26 00:08:55  wasistho
! added multiple discr.types and order
!
! Revision 1.7  2003/10/24 03:55:56  wasistho
! screen output RaNS dissipation coef k2
!
! Revision 1.6  2003/10/17 20:24:11  wasistho
! added IF statement for rans model constants
!
! Revision 1.5  2003/10/07 20:34:13  wasistho
! added screen print RaNS/DES numeric
!
! Revision 1.3  2003/08/06 15:57:43  wasistho
! added vorticities computation
!
! Revision 1.2  2003/05/31 01:47:00  wasistho
! installed turb. wall layer model
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************







