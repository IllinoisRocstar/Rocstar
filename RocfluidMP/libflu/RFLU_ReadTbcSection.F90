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
! ******************************************************************************
!
! Purpose: Read in user input related to time-dependent boundary conditions.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region type
!   tbcType        tbc type
!
! Output: None.
!
! Notes:
!   1. Sample input sections for TBCs are given in the files updateTbc*.F90
!
! ******************************************************************************
!
! $Id: RFLU_ReadTbcSection.F90,v 1.12 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadTbcSection(pRegion,tbcType)

  USE ModDataTypes
  USE ModBndPatch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY: ReadPatchSection

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: tbcType
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX = 10

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX),bcvar,pwisestr
  CHARACTER(32) :: bctitle
  CHARACTER(256) :: fname, line
  CHARACTER(CHRLEN) :: bcName,RCSIdentString
  INTEGER :: bcCode,distrib,errorFlag,ftype,iPatch,iReg,n,nbvals,nData, &
             njumps,nparams,nsvals,nswitches,nVals,prbeg,prend,var
  INTEGER :: switches(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: holdvals,morevals
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_bcvalues), POINTER :: bc
  TYPE(t_tbcvalues), POINTER :: tbc
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadTbcSection.F90,v $ $Revision: 1.12 $'

  global => pRegion%global

  CALL RegisterFunction( global,'RFLU_ReadTbcSection',&
  'RFLU_ReadTbcSection.F90' )

  IF ( global%flowType == FLOW_STEADY ) THEN
    CALL ErrorStop( global,ERR_STEADY,__LINE__ )
  END IF ! global%flowType

! ******************************************************************************
! Read BC type and variable name
! ******************************************************************************

  DO ! find first non-blank line of section
    READ(IF_INPUT,'(A256)',ERR=10,END=10) line
    IF ( LEN_TRIM(line) > 0 ) THEN
      EXIT
    END IF ! LEN_TRIM
  END DO ! <empty>

  READ(line,*) bctitle,bcvar

! ******************************************************************************
! Find the integers describing module type and variable name
! ******************************************************************************

  bcCode = -1
  ftype  = -1

  SELECT CASE( TRIM(bctitle) )

! ==============================================================================
!   Mixture BCs: ftype = FTYPE_MIXT
! ==============================================================================

! ------------------------------------------------------------------------------
!   Slip wall
! ------------------------------------------------------------------------------

    CASE ('SLIPW')
      ftype = FTYPE_MIXT
      bcCode = BC_SLIPWALL
      CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'SLIPW')

! ------------------------------------------------------------------------------
!   No-slip wall
! ------------------------------------------------------------------------------

    CASE ('NOSLIP')
      ftype = FTYPE_MIXT
      bcCode = BC_NOSLIPWALL

      SELECT CASE(TRIM(bcvar))
        CASE ('TWALL')
          var = BCDAT_NOSLIP_TWALL
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'NOSLIP')
      END SELECT

! ------------------------------------------------------------------------------
!   Inflow (based on total quantities and flow angles)
! ------------------------------------------------------------------------------

    CASE ('INFLOW_TOTANG')
      ftype = FTYPE_MIXT
      bcCode = BC_INFLOW_TOTANG

      SELECT CASE(TRIM(bcvar))
        CASE ('PTOT')
          var = BCDAT_INFLOW_PTOT
        CASE ('TTOT')
          var = BCDAT_INFLOW_TTOT
        CASE ('BETAH')
          var = BCDAT_INFLOW_BETAH
        CASE ('BETAV')
          var = BCDAT_INFLOW_BETAV
        CASE ('MACH')
          var = BCDAT_INFLOW_MACH
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'INFLOW_TOTANG')
      END SELECT

! ------------------------------------------------------------------------------
!   Inflow (based on velocity and temperature)
! ------------------------------------------------------------------------------

    CASE ('INFLOW_VELTEMP')
      ftype = FTYPE_MIXT
      bcCode = BC_INFLOW_VELTEMP

      SELECT CASE(TRIM(bcvar))
        CASE ('VELX')
          var = BCDAT_INFLOW_U
        CASE ('VELY')
          var = BCDAT_INFLOW_V
        CASE ('VELZ')
          var = BCDAT_INFLOW_W
        CASE ('TEMP')
          var = BCDAT_INFLOW_T
        CASE ('PRESS')
          var = BCDAT_INFLOW_P
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'INFLOW_VELTEMP')
      END SELECT

! ------------------------------------------------------------------------------
!   Outflow
! ------------------------------------------------------------------------------

    CASE ('OUTFLOW')
      ftype = FTYPE_MIXT
      bcCode = BC_OUTFLOW

      SELECT CASE(TRIM(bcvar))
        CASE ('PRESS')
          var = BCDAT_OUTFLOW_PRESS
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'OUTFLOW')
      END SELECT

! ------------------------------------------------------------------------------
!   Farfield
! ------------------------------------------------------------------------------

    CASE ('FARF')
      ftype = FTYPE_MIXT
      bcCode = BC_FARFIELD

      SELECT CASE(TRIM(bcvar))
        CASE ('MACH')
          var = BCDAT_FARF_MACH
        CASE ('ATTACK')
          var = BCDAT_FARF_ATTACK
        CASE ('SLIP')
          var = BCDAT_FARF_SLIP
        CASE ('PRESS')
          var = BCDAT_FARF_PRESS
        CASE ('TEMP')
          var = BCDAT_FARF_TEMP
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'FARF')
      END SELECT

! ------------------------------------------------------------------------------
!   Injection
! ------------------------------------------------------------------------------

    CASE ('INJECT')
      ftype = FTYPE_MIXT
      bcCode = BC_INJECTION

      SELECT CASE(TRIM(bcvar))
        CASE ('MFRATE')
          var = BCDAT_INJECT_MFRATE
        CASE ('TEMP')
          var = BCDAT_INJECT_TEMP
        CASE DEFAULT
          CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'INJECT')
      END SELECT

#ifdef TURB
! ==============================================================================
!   Turbulence BCs: ftype = FTYPE_TURB
! ==============================================================================

! TO DO
! END TO DO
#endif

#ifdef SPEC
! ==============================================================================
!   Species BCs: ftype = FTYPE_SPEC
! ==============================================================================

! ------------------------------------------------------------------------------
!   Inflow
! ------------------------------------------------------------------------------

    CASE ('SPEC_INFLOW')
      ftype = FTYPE_SPEC
      bcCode = BC_INFLOW

      IF ( bcvar(1:4) == 'SPEC' ) THEN
        IF ( LEN_TRIM(bcvar) == 4 ) THEN
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '*** WARNING *** Improper specification in SPEC_INFLOW TBC.'
          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '                Assuming SPEC to mean SPEC1.'

          var = 1
        ELSE
          READ(bcvar(5:),*) var
        END IF ! LEN_TRIM
      ELSE
        CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'SPEC_INFLOW')
      END IF ! bcvar

! ------------------------------------------------------------------------------
!   Farfield
! ------------------------------------------------------------------------------

    CASE ('SPEC_FARF')
      ftype = FTYPE_SPEC
      bcCode = BC_FARFIELD

      IF ( bcvar(1:4) == 'SPEC' ) THEN
        IF ( LEN_TRIM(bcvar) == 4 ) THEN
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '*** WARNING *** Improper specification in SPEC_FARF TBC.'
          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '                Assuming SPEC to mean SPEC1.'

          var = 1
        ELSE
          READ(bcvar(5:),*) var
        END IF ! LEN_TRIM
      ELSE
        CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'SPEC_FARF')
      END IF ! bcvar

! ------------------------------------------------------------------------------
!   Injection
! ------------------------------------------------------------------------------

    CASE ('SPEC_INJECT')
      ftype = FTYPE_SPEC
      bcCode = BC_INJECTION

      IF ( bcvar(1:4) == 'SPEC' ) THEN
        IF ( LEN_TRIM(bcvar) == 4 ) THEN
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '*** WARNING *** Improper specification in SPEC_INJECT TBC.'
          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                '                Assuming SPEC to mean SPEC1.'

          var = 1
        ELSE
          READ(bcvar(5:),*) var
        END IF ! LEN_TRIM
      ELSE
        CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__,'SPEC_INJECT')
      END IF ! bcvar
#endif

#ifdef RADI
! ==============================================================================
!   Radiation BCs: ftype = FTYPE_RADI
! ==============================================================================

! TO DO
! END TO DO
#endif

    CASE DEFAULT
      CALL ErrorStop(global,ERR_UNKNOWN_BC,__LINE__)
  END SELECT ! TRIM(bctitle)

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

  keys(TBCDAT_ONTIME)  = 'ONTIME'
  keys(TBCDAT_OFFTIME) = 'OFFTIME'
  keys(TBCDAT_AMP)     = 'AMP'

  SELECT CASE ( tbcType )
    CASE ( TBC_SINUSOIDAL )
      keys(TBCDAT_FREQ)  = 'FREQ'
      keys(TBCDAT_PHASE) = 'PHASE'
      nVals = 5
    CASE ( TBC_STOCHASTIC )
      keys(TBCDAT_TIMECOR) = 'TIMECOR'
      keys(TBCDAT_SHAPE)   = 'SHAPE'
      keys(TBCDAT_MINCUT)  = 'MINCUT'
      keys(TBCDAT_MAXCUT)  = 'MAXCUT'
      nVals = 7
    CASE ( TBC_WHITENOISE )
      keys(4) = 'SUBSTEP' ! becomes a switch
      nVals = 4
    CASE ( TBC_PIECEWISE )
      keys(3) = 'ORDER'   ! becomes a switch (overwrites 'AMP')
      keys(4) = 'NJUMPS'  ! becomes a switch
      nVals = 4
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! tbcType

  CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,prbeg,prend, &
                        distrib,fname,bcName,defined)

  IF (.NOT. defined(TBCDAT_ONTIME))      vals(TBCDAT_ONTIME)  = -1.E20_RFREAL
  IF (vals(TBCDAT_ONTIME) < 0.0_RFREAL)  vals(TBCDAT_ONTIME)  = -1.E20_RFREAL

  IF (.NOT. defined(TBCDAT_OFFTIME))     vals(TBCDAT_OFFTIME) =  1.E20_RFREAL
  IF (vals(TBCDAT_OFFTIME) < 0.0_RFREAL) vals(TBCDAT_OFFTIME) =  1.E20_RFREAL

  IF (tbcType /= TBC_PIECEWISE) THEN
    IF (.NOT. defined(TBCDAT_AMP)) &
      CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)

    IF (vals(TBCDAT_AMP) <= 0.0_RFREAL) &
      CALL ErrorStop(global,ERR_VAL_BCVAL,__LINE__)
  END IF ! tbcType

! ******************************************************************************
! Set values of vals and switches
! ******************************************************************************

  SELECT CASE ( tbcType )

! ==============================================================================
!   Sinusoidal variation
! ==============================================================================

    CASE ( TBC_SINUSOIDAL )
      nswitches = 0
      nparams = 5
      nsvals = 1
      nbvals = 0

      IF (.NOT. defined(TBCDAT_FREQ)) &
        CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)

      IF (.NOT. defined(TBCDAT_PHASE)) vals(TBCDAT_PHASE) = 0.0_RFREAL

      vals(TBCDAT_PHASE) = vals(TBCDAT_PHASE)*global%deg2rad

! ==============================================================================
!   Stochastic variation
! ==============================================================================

    CASE ( TBC_STOCHASTIC )
      nswitches = 0
      nparams = 7
      nsvals = 0
      nbvals = 3

      IF (.NOT. defined(TBCDAT_TIMECOR)) &
        CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)

      IF (vals(TBCDAT_TIMECOR) <= 0.0_RFREAL) &
        CALL ErrorStop(global,ERR_VAL_BCVAL,__LINE__)

      IF (.NOT. defined(TBCDAT_SHAPE)) vals(TBCDAT_SHAPE) = 1.0_RFREAL

      IF (vals(TBCDAT_SHAPE) < 0.1_RFREAL .OR. vals(TBCDAT_SHAPE) > 10.0_RFREAL) &
        CALL ErrorStop(global,ERR_VAL_BCVAL,__LINE__)

      IF (.NOT. defined(TBCDAT_MINCUT)) vals(TBCDAT_MINCUT) = -1.0_RFREAL

      IF (.NOT. defined(TBCDAT_MAXCUT)) vals(TBCDAT_MAXCUT) = -1.0_RFREAL

! ==============================================================================
!   Whitenoise variation
! ==============================================================================

    CASE ( TBC_WHITENOISE )
      nswitches = 1
      nparams = 3
      nsvals = 0
      nbvals = 1

      IF (.NOT. defined(4)) switches(TBCSWI_SUBSTEP) = TBCOPT_STEP

      IF ( vals(4) < 0.1_RFREAL ) THEN
        switches(TBCSWI_SUBSTEP) = TBCOPT_STEP
      ELSE
        switches(TBCSWI_SUBSTEP) = TBCOPT_SUBSTEP
      END IF ! vals

! ==============================================================================
!   Piecewise variation
! ==============================================================================

    CASE ( TBC_PIECEWISE )
      nswitches = 2
      nsvals = 1
      nbvals = 0

      IF (.NOT. defined(3)) switches(TBCSWI_ORDER) = TBCOPT_CONSTANT

      IF (vals(3) < 0.1_RFREAL) THEN
        switches(TBCSWI_ORDER) = TBCOPT_CONSTANT
      ELSE
        switches(TBCSWI_ORDER) = TBCOPT_LINEAR
      END IF ! vals

      IF (.NOT. defined(4)) THEN
        CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
      END IF

      njumps = NINT(vals(4))
      switches(TBCSWI_NJUMPS) = njumps

      IF (njumps < 1 .OR. njumps > 1000000) &
        CALL ErrorStop(global,ERR_VAL_BCVAL,__LINE__)

      SELECT CASE (switches(TBCSWI_ORDER))
        CASE (TBCOPT_CONSTANT)
          nparams = TBCDAT_DAT0 + 2*njumps+1
        CASE (TBCOPT_LINEAR)
          nparams = TBCDAT_DAT0 + 2*njumps
      END SELECT ! switches(TBCSWI_ORDER)

      ALLOCATE(morevals(nparams),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

      morevals(1:TBCDAT_DAT0) = vals(1:TBCDAT_DAT0)

      DO n=TBCDAT_DAT0+1,nparams
        DO ! find first non-blank line
          READ(IF_INPUT,'(A256)',ERR=10,END=10) line
          IF (LEN_TRIM(line) > 0) EXIT
        END DO ! <empty>

        IF (line(1:1) == '#') GOTO 10

        READ(line,*) pwisestr, morevals(n)
      END DO ! n

! - check that the times entered are increasing

      DO n=TBCDAT_DAT0+4,nparams,2
        IF (morevals(n) .le. morevals(n-2)) &
          CALL ErrorStop(global,ERR_VAL_BCVAL,__LINE__)
      END DO ! n

      DO ! find end of section
        READ(IF_INPUT,'(A256)',ERR=10,END=10) line
        IF (line(1:1) == '#') EXIT
      END DO ! <empty>

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! tbctype

! ******************************************************************************
! Copy values to variables
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================
!   Check whether this global patch exists in this region
! ==============================================================================

    IF ( (pPatch%bcType >= bcCode .AND. &
          pPatch%bcType <= bcCode+BC_RANGE) .AND. &
         (pPatch%iPatchGlobal >= prbeg .AND. &
          pPatch%iPatchGlobal <= prend) ) THEN
      SELECT CASE (ftype)
        CASE (FTYPE_MIXT)
          bc => pPatch%mixt
        CASE (FTYPE_TURB)
          bc => pPatch%turb
        CASE (FTYPE_SPEC)
          bc => pPatch%spec
        CASE (FTYPE_RADI)
          bc => pPatch%valRadi
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ftype

      nData =  bc%nData
      tbc   => bc%tbcs(var)

      IF ( var > nData ) THEN
        CALL ErrorStop(global,ERR_BC_VARNAME,__LINE__)
      END IF ! var

      IF ( tbc%tbcType /= TBC_NONE ) THEN
        CALL ErrorStop(global,ERR_PATCH_OVERSPEC,__LINE__,'TBC')
      END IF ! tbc%tbcType

      tbc%tbcType = tbcType

      IF ( nswitches > 0 ) THEN
        ALLOCATE(tbc%switches(nswitches),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        tbc%switches = switches(1:nswitches)
      END IF ! nswitches

      IF ( nparams > 0 ) THEN
        ALLOCATE(tbc%params(nparams),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        IF (tbcType == TBC_PIECEWISE) THEN
          tbc%params = morevals(1:nparams)
        ELSE
          tbc%params = vals(1:nparams)
        END IF ! tbcType
      END IF ! nparams

      IF (nsvals > 0) THEN
        ALLOCATE(tbc%svals(nsvals),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        tbc%svals = 0.0_RFREAL
      END IF ! nsvals

      IF (nbvals > 0) THEN
        ALLOCATE(tbc%bvals(nbvals,pPatch%nBFaces),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        tbc%bvals = 0.0_RFREAL

        IF ( tbcType == TBC_STOCHASTIC ) THEN
          tbc%bvals(TBCSTO_VAL,:) = 0.5_RFREAL * tbc%params(TBCDAT_AMP)**2
        END IF ! tbcType
      END IF ! nbvals

      distrib = bc%distrib ! original value of bc%distrib

      IF (tbcType == TBC_STOCHASTIC .OR. & ! tbc types forcing BCDAT_DISTRIB
        tbcType == TBC_WHITENOISE) bc%distrib = BCDAT_DISTRIB

      IF ( bc%distrib == BCDAT_DISTRIB ) THEN
        ALLOCATE(tbc%mean(pPatch%nBFaces),STAT=errorFlag)
      ELSE
        ALLOCATE(tbc%mean(0:1),STAT=errorFlag)
      END IF ! bc%distrib

      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

! ==============================================================================
!     If distrib needs to change for BC, allocate and fill expanded bc%vals
!     arrays
! ==============================================================================

      IF ( bc%distrib == BCDAT_DISTRIB .AND. distrib == BCDAT_CONSTANT ) THEN
        ALLOCATE( holdvals(nData),STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        holdvals = bc%vals(:,1)

        DEALLOCATE(bc%vals,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
        END IF ! global%error

        ALLOCATE(bc%vals(nData,pPatch%nBFaces),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global%error

        DO n = 1, nData
          bc%vals(n,:) = holdvals(n)
        END DO ! n

        DEALLOCATE (holdvals,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
        END IF ! global%error

      END IF ! bc%distrib

      tbc%mean = bc%vals(var,:)

      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

    END IF ! patch in region, correct boundary type
  END DO ! iPatch

  IF ( ALLOCATED(morevals) .EQV. .TRUE. ) THEN
    DEALLOCATE(morevals,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
    END IF ! global%error
  END IF ! ALLOCATED

  GOTO 999

! ******************************************************************************
! Error handling
! ******************************************************************************

10  CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__)

999 CONTINUE

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_ReadTbcSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadTbcSection.F90,v $
! Revision 1.12  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/08/19 15:38:48  mparmar
! Renamed patch variables
!
! Revision 1.9  2005/04/27 02:09:00  haselbac
! Added support for INFLOW_VELTEMP
!
! Revision 1.8  2004/10/19 19:25:17  haselbac
! Removed RFVFx on injecting boundaries, clean-up
!
! Revision 1.7  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.6  2004/02/26 21:01:54  haselbac
! Removed stray TEMPORARY statement
!
! Revision 1.5  2004/01/29 22:56:46  haselbac
! Added species stuff, clean-up
!
! Revision 1.4  2003/07/22 01:59:42  haselbac
! Added global%warnCounter
!
! Revision 1.3  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.2  2003/06/04 20:06:38  jferry
! added capability for PEUL BCs with TBCs
!
! Revision 1.1  2002/11/16 00:05:31  haselbac
! Moved here from rocflu directory
!
! Revision 1.4  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.3  2002/09/28 14:33:24  jferry
! added relative velocity to injection BC
!
! Revision 1.2  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.1  2002/09/25 18:31:57  jferry
! Implemented TBCs for ROCFLU
!
! ******************************************************************************







