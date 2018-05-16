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
! Purpose: Set values derived from user input.
!
! Description: Derived variables are mostly components of turbInput data type.
!
! Input: regions = Input parameters for all regions.
!
! Output: regions = Derived variables stored as turbInput data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_DerivedInputValues.F90,v 1.11 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_DerivedInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModTurbulence, ONLY : t_turb_input
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, m

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER     :: global
  TYPE(t_turb_input), POINTER :: input

  INTEGER :: flowModel, turbModel, iLev, nSv, nSt, errorFlag
#ifdef RFLO
  REAL(RFREAL)      :: one3rd, filtLScale, rNdel(3)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_DerivedInputValues.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_DerivedInputValues',&
  'TURB_DerivedInputValues.F90' )

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Entering TURB_DerivedInputValues...'
  END IF ! global%verbLevel

! set local constants ---------------------------------------------------------

#ifdef RFLO
  one3rd = 1._RFREAL/3._RFREAL
#endif

! region related data (all levels) --------------------------------------------

#ifdef RFLO
  DO iReg = 1,global%nRegions
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif

    flowModel =  regions(iReg)%mixtInput%flowModel
    turbModel =  regions(iReg)%mixtInput%turbModel
    input     => regions(iReg)%turbInput

! - number of N-S transport variables

#ifdef RFLO
    DO iLev=1,regions(iReg)%nGridLevels
      regions(iReg)%levels(iLev)%mixt%nTv = 4 
    ENDDO
#endif
#ifdef RFLU
    regions(iReg)%mixtInput%nTv = 4 
#endif

! - number of dummy cells

#ifdef RFLO
    IF (regions(iReg)%nDumCells < 3) THEN
      IF ((turbModel==TURB_MODEL_SCALSIM).OR. &
          (turbModel==TURB_MODEL_DYNSMAG)) THEN     ! LES with filtering
        IF ((input%filterWidth(DIRI) > 1) .OR. &
            (input%filterWidth(DIRJ) > 1) .OR. &
            (input%filterWidth(DIRK) > 1)) THEN
          regions(iReg)%nDumCells = 3
        ENDIF
      ELSEIF (turbModel==TURB_MODEL_DYNMIXD .OR. &
              turbModel==TURB_MODEL_HDESSA) THEN
        regions(iReg)%nDumCells = 3
      ELSE
      ENDIF
    ENDIF
#endif

! - GENERAL TURB --------------------------------------------------------------
! - class of turbulence model

    IF ((turbModel==TURB_MODEL_FIXSMAG).OR. &
        (turbModel==TURB_MODEL_SCALSIM).OR. &
        (turbModel==TURB_MODEL_DYNSMAG).OR. &
        (turbModel==TURB_MODEL_DYNMIXD)) THEN   
      input%modelClass = MODEL_LES
    ELSEIF ((turbModel==TURB_MODEL_SA).OR. &
            (turbModel==TURB_MODEL_DESSA).OR. &                 
            (turbModel==TURB_MODEL_HDESSA)) THEN                 
      input%modelClass = MODEL_RANS
    ENDIF

! - number of conservative and derived variables

    input%nCv = 0
    input%nDv = 0
    IF (input%modelClass == MODEL_LES) THEN
      input%nDv = 1
    ELSEIF (input%modelClass == MODEL_RANS) THEN
      IF ((turbModel==TURB_MODEL_SA).OR. & 
          (turbModel==TURB_MODEL_DESSA).OR. &
          (turbModel==TURB_MODEL_HDESSA)) THEN
        input%nCv = 1
      ENDIF
    ENDIF

! - number of zero-one switch fields

    IF (turbModel==TURB_MODEL_FIXSMAG .OR. &
        turbModel==TURB_MODEL_DYNSMAG .OR. &
        turbModel==TURB_MODEL_DYNMIXD) THEN
      input%nZof = input%nZof + 1
      global%calcFaceCtr = .TRUE.
    ENDIF

! - number of permanent and collected variables to be time averaged

    input%nFixSt = 3
    nSt = 0
#ifdef STATS
    IF ((global%flowType == FLOW_UNSTEADY) .AND. &
        (global%doStat == ACTIVE) .AND. &
        (global%turbNStat > 0)) THEN
      nSt = 3
    ENDIF
#endif
    input%nSt = 0  ! = nSt to activate

! - number of stress components

    nSv       = 6
    input%nSv = 0  ! = nSv to activate

! - number of gradient variables

    IF ((turbModel == TURB_MODEL_DYNSMAG).OR. &
        (turbModel == TURB_MODEL_DYNMIXD)) THEN
#ifdef RFLO
      input%nGrad = 9
#endif
#ifdef RFLU
      input%nGrad = 3
#endif
    ELSEIF ((turbModel == TURB_MODEL_SA).OR. &
            (turbModel == TURB_MODEL_DESSA).OR. &
            (turbModel == TURB_MODEL_HDESSA)) THEN
#ifdef RFLO
      input%nGrad = 3
#endif
#ifdef RFLU
      input%nGrad = 1
#endif
    ELSE
      input%nGrad = 0
    ENDIF

! - LES SPECIFIC --------------------------------------------------------------
! - filter width scaling factor

#ifdef RFLO
    IF (input%deltaType == DELTYPE_CBRT) THEN
      filtLScale = 1._RFREAL
    ELSEIF (input%deltaType == DELTYPE_SQRT) THEN
      filtLScale = 3._RFREAL
    ENDIF

    rNdel(1) = REAL(input%filterWidth(DIRI)) 
    rNdel(2) = REAL(input%filterWidth(DIRJ)) 
    rNdel(3) = REAL(input%filterWidth(DIRK)) 
    input%delFac2 = (rNdel(1)*rNdel(2)*rNdel(3))**one3rd
    IF (input%delFac2 < REAL_SMALL) THEN
      input%delFac2 = SQRT(rNdel(1)*rNdel(2))+ &
                      SQRT(rNdel(1)*rNdel(3))+ &
                      SQRT(rNdel(2)*rNdel(3))
    ENDIF
    IF ((input%filterWidth(DIRI)/=0.AND.input%filterWidth(DIRJ)==0      &
                                   .AND.input%filterWidth(DIRK)==0).OR. &
        (input%filterWidth(DIRI)==0.AND.input%filterWidth(DIRJ)/=0      &
                                   .AND.input%filterWidth(DIRK)==0).OR. &
        (input%filterWidth(DIRI)==0.AND.input%filterWidth(DIRJ)==0      &
                                   .AND.input%filterWidth(DIRK)/=0)) THEN
      input%delFac2 = rNdel(1)+rNdel(2)+rNdel(3)
    ENDIF
    input%delFac2 = filtLScale*input%delFac2**2
#endif
#ifdef RFLU
    input%delFac2 = REAL(input%filterWidth(DIRI))
    input%delFac2 = input%delFac2**2
#endif

    IF (turbModel==TURB_MODEL_FIXSMAG) THEN
      input%delFac2 = 1._RFREAL
    ENDIF

! - RANS SPECIFIC -------------------------------------------------------------
! - model constants

#ifdef RFLO
    IF (global%flowType == FLOW_UNSTEADY) THEN
      input%smoocf = -1._RFREAL
    ENDIF
#endif

    IF ((turbModel == TURB_MODEL_SA).OR. &
        (turbModel == TURB_MODEL_DESSA).OR. &
        (turbModel == TURB_MODEL_HDESSA)) THEN
      ALLOCATE( input%const(MC_SA_NELM),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      input%const(MC_SA_CB1)   = 0.1355_RFREAL
      input%const(MC_SA_CB2)   = 0.622_RFREAL
      input%const(MC_SA_CW2)   = 0.3_RFREAL
      input%const(MC_SA_CW3)   = 2.0_RFREAL
      input%const(MC_SA_CV1)   = 7.1_RFREAL
      input%const(MC_SA_RSIG)  = 3.0_RFREAL/2.0_RFREAL
      input%const(MC_SA_RKAP)  = 1.0_RFREAL/0.41_RFREAL
      input%const(MC_SA_CW1)   = &
            input%const(MC_SA_CB1)*input%const(MC_SA_RKAP)**2 + &
            (1._RFREAL+input%const(MC_SA_CB2))*input%const(MC_SA_RSIG)
    ELSE
      NULLIFY( input%const )
    ENDIF

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Leaving TURB_DerivedInputValues.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_DerivedInputValues.F90,v $
! Revision 1.11  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/02/04 05:00:05  wasistho
! added enter and leave statements
!
! Revision 1.8  2006/01/17 17:25:22  wasistho
! applied tripping to all eddy viscosity models
!
! Revision 1.7  2006/01/12 09:49:56  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.6  2005/04/18 21:40:02  wasistho
! use nDumCell=3 for hybrid DES
!
! Revision 1.5  2005/03/09 06:34:51  wasistho
! incorporated HDESSA
!
! Revision 1.4  2004/10/22 23:16:35  wasistho
! set max. number of collected extra statistics, nSt to 3
!
! Revision 1.3  2004/08/09 18:38:40  wasistho
! fixed input%delFac2 for RFLU
!
! Revision 1.2  2004/03/19 02:45:21  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.10  2003/10/26 00:24:52  wasistho
! another slightly more accurate model constant
!
! Revision 1.9  2003/10/24 20:56:49  wasistho
! corrected SA model constant cb2
!
! Revision 1.8  2003/10/07 23:56:05  wasistho
! set RaNS smoocf to -1.0 for unsteady flow
!
! Revision 1.6  2003/08/01 22:17:25  wasistho
! prepared rocturb for Genx
!
! Revision 1.5  2003/07/22 02:59:09  wasistho
! prepare more accurate rocturb restart
!
! Revision 1.4  2003/05/24 02:41:56  wasistho
! sv collector turned off
!
! Revision 1.3  2003/05/24 02:10:51  wasistho
! turbulence statistics expanded
!
! Revision 1.2  2002/10/15 00:00:00  wasistho
! Removed temporary safety
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







