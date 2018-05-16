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
! Purpose: read in user input related to time-dependent boundary conditions.
!
! Description: none.
!
! Input: type of TBC boundary condtion
!
! Output: regions will be updated with TBC data
!
! Notes:
!
!   Sample input sections for TBCs are given in the files updateTbc*.F90
!
!******************************************************************************
!
! $Id: RFLO_ReadTbcSection.F90,v 1.6 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadTbcSection( regions,tbcType )

  USE ModDataTypes
  USE ModBndPatch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
#ifdef PEUL
  USE PEUL_ModParameters
#endif
  USE ModInterfaces, ONLY : ReadPatchSection
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER, INTENT(IN)     :: tbcType

! ... loop variables
  INTEGER :: iReg, iPatch, n

! ... local variables
  INTEGER, PARAMETER :: NVALS_MAX = 10

  CHARACTER(10)  :: keys(NVALS_MAX), bcvar, pwisestr
  CHARACTER(32)  :: bctitle
  CHARACTER(256) :: fname, line

  INTEGER :: nVals, brbeg, brend, prbeg, prend, distrib
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag
  INTEGER :: ftype, bcCode, var, nswitches, nparams, nsvals, nbvals, nData
  INTEGER :: switches(NVALS_MAX), njumps

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL)              :: vals(NVALS_MAX)
  REAL(RFREAL), ALLOCATABLE :: holdvals(:), morevals(:)

  TYPE(t_patch),     POINTER :: patch
  TYPE(t_bcvalues),  POINTER :: bc
  TYPE(t_tbcvalues), POINTER :: tbc
  TYPE(t_global),    POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadTbcSection',&
  'RFLO_ReadTbcSection.F90' )

  IF (global%flowType == FLOW_STEADY) &
    CALL ErrorStop( global,ERR_STEADY,__LINE__ )

! read BC type and variable name ----------------------------------------------

  DO ! find first non-blank line of section
    READ(IF_INPUT,'(A256)',err=10,end=10) line
    IF (LEN_TRIM(line) > 0) EXIT
  ENDDO

  READ(line,*) bctitle, bcvar

! find the integers describing module type and variable name

  bcCode = -1
  ftype  = -1

  SELECT CASE(TRIM(bctitle))

! *** Mixture BCs *** : ftype = FTYPE_MIXT

  CASE ('SLIPW')
    ftype = FTYPE_MIXT
    bcCode = BC_SLIPWALL
    CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'SLIPW' )

  CASE ('NOSLIP')
    ftype = FTYPE_MIXT
    bcCode = BC_NOSLIPWALL

    SELECT CASE(TRIM(bcvar))
    CASE ('TWALL')
      var = BCDAT_NOSLIP_TWALL
    CASE DEFAULT
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'NOSLIP' )
    END SELECT

  CASE ('INFLOW_TOTANG')
    ftype = FTYPE_MIXT
    bcCode = BC_INFLOW

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
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'INFLOW' )
    END SELECT

  CASE ('INFLOW_VELTEMP')
    ftype = FTYPE_MIXT
    bcCode = BC_INFLOW

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
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'INFLOW_VELTEMP' )
    END SELECT

  CASE ('INFLOW_VELPRESS')
    ftype = FTYPE_MIXT
    bcCode = BC_INFLOW

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
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'INFLOW_VELPRESS' )
    END SELECT

  CASE ('OUTFLOW')
    ftype = FTYPE_MIXT
    bcCode = BC_OUTFLOW

    SELECT CASE(TRIM(bcvar))
    CASE ('PRESS')
      var = BCDAT_OUTFLOW_PRESS
    CASE DEFAULT
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'OUTFLOW' )
    END SELECT

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
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'FARF' )
    END SELECT

  CASE ('INJECT')
    ftype = FTYPE_MIXT
    bcCode = BC_INJECTION

    SELECT CASE(TRIM(bcvar))
    CASE ('MFRATE')
      var = BCDAT_INJECT_MFRATE
    CASE ('TEMP')
      var = BCDAT_INJECT_TEMP
    CASE ('RFVFU')
      var = BCDAT_INJECT_RFVFU
    CASE ('RFVFV')
      var = BCDAT_INJECT_RFVFV
    CASE ('RFVFW')
      var = BCDAT_INJECT_RFVFW
    CASE DEFAULT
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'INJECT' )
    END SELECT

! *** Turbulence BCs *** : ftype = FTYPE_TURB
! *** Species BCs *** : ftype = FTYPE_SPEC

! *** Eulerian Particle BCs *** : ftype = FTYPE_PEUL

#ifdef PEUL
  CASE ('PEUL_INFLOW')
    ftype = FTYPE_PEUL
    bcCode = BC_INFLOW

    IF (bcvar(1:4) == 'DENS') THEN
      IF (LEN_TRIM(bcvar) == 4) THEN
        WRITE(*,*) '### WARNING: Improper specification in PEUL_INFLOW TBC'
        WRITE(*,*) '             assuming DENS to mean DENS1'
        var = 1
      ELSE
        READ(bcvar(5:),*) var
      END IF
    ELSE
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'PEUL_INFLOW' )
    END IF

  CASE ('PEUL_FARF')
    ftype = FTYPE_PEUL
    bcCode = BC_FARFIELD

    IF (bcvar(1:4) == 'DENS') THEN
      IF (LEN_TRIM(bcvar) == 4) THEN
        WRITE(*,*) '### WARNING: Improper specification in PEUL_FARF TBC'
        WRITE(*,*) '             assuming DENS to mean DENS1'
        var = 1
      ELSE
        READ(bcvar(5:),*) var
      END IF
    ELSE
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'PEUL_FARF' )
    END IF

  CASE ('PEUL_INJECT')
    ftype = FTYPE_PEUL
    bcCode = BC_INJECTION

    IF (bcvar(1:4) == 'FRAC') THEN
      IF (LEN_TRIM(bcvar) == 4) THEN
        WRITE(*,*) '### WARNING: Improper specification in PEUL_INJECT TBC'
        WRITE(*,*) '             assuming FRAC to mean FRAC1'
        var = 1
      ELSE
        READ(bcvar(5:),*) var
      END IF
    ELSE
      CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__,'PEUL_INJECT' )
    END IF
#endif

! *** Radiation BCs *** : ftype = FTYPE_RADI

  CASE DEFAULT
    CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )

  END SELECT

! specify keywords and search for them ----------------------------------------

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
    CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT ! tbcType

  CALL ReadPatchSection( global,IF_INPUT,nVals,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )

  IF (.NOT. defined(TBCDAT_ONTIME))      vals(TBCDAT_ONTIME)  = -1.E20_RFREAL
  IF (vals(TBCDAT_ONTIME) < 0.0_RFREAL)  vals(TBCDAT_ONTIME)  = -1.E20_RFREAL

  IF (.NOT. defined(TBCDAT_OFFTIME))     vals(TBCDAT_OFFTIME) =  1.E20_RFREAL
  IF (vals(TBCDAT_OFFTIME) < 0.0_RFREAL) vals(TBCDAT_OFFTIME) =  1.E20_RFREAL

  IF (tbcType /= TBC_PIECEWISE) THEN
    IF (.NOT. defined(TBCDAT_AMP)) &
      CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

    IF (vals(TBCDAT_AMP) <= 0.0_RFREAL) &
      CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__ )
  ENDIF ! tbcType

! set values of vals and switches ---------------------------------------------

  SELECT CASE ( tbcType )

  CASE ( TBC_SINUSOIDAL )

    nswitches = 0
    nparams   = 5
    nsvals    = 1
    nbvals    = 0

    IF (.NOT. defined(TBCDAT_FREQ)) &
      CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

    IF (.NOT. defined(TBCDAT_PHASE)) vals(TBCDAT_PHASE) = 0.0_RFREAL
    vals(TBCDAT_PHASE) = vals(TBCDAT_PHASE)*global%rad

  CASE ( TBC_STOCHASTIC )

    nswitches = 0
    nparams   = 7
    nsvals    = 0
    nbvals    = 3

    IF (.NOT. defined(TBCDAT_TIMECOR)) &
      CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

    IF (vals(TBCDAT_TIMECOR) <= 0.0_RFREAL) &
      CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__ )

    IF (.NOT. defined(TBCDAT_SHAPE)) vals(TBCDAT_SHAPE) = 1.0_RFREAL

    IF (vals(TBCDAT_SHAPE) < 0.1_RFREAL .OR. vals(TBCDAT_SHAPE) > 10.0_RFREAL) &
      CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__ )

    IF (.NOT. defined(TBCDAT_MINCUT)) vals(TBCDAT_MINCUT) = -1.0_RFREAL
    IF (.NOT. defined(TBCDAT_MAXCUT)) vals(TBCDAT_MAXCUT) = -1.0_RFREAL

  CASE ( TBC_WHITENOISE )

    nswitches = 1
    nparams = 3
    nsvals = 0
    nbvals = 1

    IF (.NOT. defined(4)) switches(TBCSWI_SUBSTEP) = TBCOPT_STEP
    IF (vals(4) < 0.1_RFREAL) THEN
      switches(TBCSWI_SUBSTEP) = TBCOPT_STEP
    ELSE
      switches(TBCSWI_SUBSTEP) = TBCOPT_SUBSTEP
    ENDIF

  CASE ( TBC_PIECEWISE )

    nswitches = 2
    nsvals = 1
    nbvals = 0

    IF (.NOT. defined(3)) switches(TBCSWI_ORDER) = TBCOPT_CONSTANT
    IF (vals(3) < 0.1_RFREAL) THEN
      switches(TBCSWI_ORDER) = TBCOPT_CONSTANT
    ELSE
      switches(TBCSWI_ORDER) = TBCOPT_LINEAR
    ENDIF

    IF (.NOT. defined(4)) &
      CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )
    njumps = NINT(vals(4))
    switches(TBCSWI_NJUMPS) = njumps
    IF (njumps < 1 .OR. njumps > 1000000) &
      CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__ )

    SELECT CASE (switches(TBCSWI_ORDER))
    CASE (TBCOPT_CONSTANT)
      nparams = TBCDAT_DAT0 + 2*njumps+1
    CASE (TBCOPT_LINEAR)
      nparams = TBCDAT_DAT0 + 2*njumps
    END SELECT ! switches(TBCSWI_ORDER)

    ALLOCATE( morevals(nparams),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    morevals(1:TBCDAT_DAT0) = vals(1:TBCDAT_DAT0)
    DO n=TBCDAT_DAT0+1,nparams
      DO ! find first non-blank line
        READ(IF_INPUT,'(A256)',err=10,end=10) line
        IF (LEN_TRIM(line) > 0) EXIT
      ENDDO
      IF (line(1:1) == '#') GOTO 10

      READ(line,*) pwisestr, morevals(n)
    ENDDO ! n

! - check that the times entered are increasing

    DO n=TBCDAT_DAT0+4,nparams,2
      IF (morevals(n) .le. morevals(n-2)) &
        CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__ )
    ENDDO ! n

    DO ! find end of section
      READ(IF_INPUT,'(A256)',err=10,end=10) line
      IF (line(1:1) == '#') EXIT
    ENDDO

  CASE DEFAULT
    CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT

! copy values to variables -----------------------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=bcCode .AND. &
           patch%bcType<=bcCode+BC_RANGE)       .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        n1    = ABS(patch%l1end-patch%l1beg)
        n2    = ABS(patch%l2end-patch%l2beg)
        iOff  = n1 + 1
        ijBeg = IndIJ( 0, 0,iOff)
        ijEnd = IndIJ(n1,n2,iOff)

        SELECT CASE (ftype)

        CASE (FTYPE_MIXT)
          bc => patch%mixt
        CASE (FTYPE_TURB)
          bc => patch%turb
        CASE (FTYPE_PEUL)
          bc => patch%peul
        CASE (FTYPE_SPEC)
          bc => patch%spec
        CASE (FTYPE_RADI)
          bc => patch%valRadi
        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

        END SELECT

        nData = bc%nData
        tbc  => bc%tbcs(var)

        IF (.NOT. bc%bcSet) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__ )

        IF (var > nData) THEN
          WRITE(*,*) 'ftype, bcCode, nData, var: ', ftype, bcCode, nData, var
          CALL ErrorStop( global,ERR_BC_VARNAME,__LINE__ )
        ENDIF

        IF (tbc%tbcType /= TBC_NONE) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'TBC' )

        tbc%tbcType = tbcType

        IF (nswitches > 0) THEN
          ALLOCATE( tbc%switches(nswitches),stat=errorFlag )
          tbc%switches = switches(1:nswitches)
        ENDIF

        IF (nparams > 0) THEN
          ALLOCATE( tbc%params(nparams),stat=errorFlag )
          IF (tbcType == TBC_PIECEWISE) THEN
            tbc%params = morevals(1:nparams)
          ELSE
            tbc%params = vals(1:nparams)
          ENDIF
        ENDIF

        IF (nsvals > 0) THEN
          ALLOCATE( tbc%svals(nsvals),stat=errorFlag )
          tbc%svals = 0.0_RFREAL
        ENDIF

        IF (nbvals > 0) THEN
          ALLOCATE( tbc%bvals(nbvals,ijBeg:ijEnd),stat=errorFlag )
          tbc%bvals = 0.0_RFREAL
          IF (tbcType == TBC_STOCHASTIC) &
            tbc%bvals(TBCSTO_VAL,:) = 0.5_RFREAL * tbc%params(TBCDAT_AMP)**2
        ENDIF

        distrib = bc%distrib ! original value of bc%distrib

        IF (tbcType == TBC_STOCHASTIC .OR. & ! tbc types forcing BCDAT_DISTRIB
            tbcType == TBC_WHITENOISE) bc%distrib = BCDAT_DISTRIB

        IF (bc%distrib == BCDAT_DISTRIB) THEN
          ALLOCATE( tbc%mean(ijBeg:ijEnd),stat=errorFlag )
        ELSE
          ALLOCATE( tbc%mean(0:1),stat=errorFlag )
        ENDIF

        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- if distrib needs to change for BC, allocate and fill
!       expanded bc%vals arrays

        IF (bc%distrib == BCDAT_DISTRIB .AND. distrib == BCDAT_CONSTANT) THEN

          ALLOCATE( holdvals(nData),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          holdvals = bc%vals(:,1)

          DEALLOCATE( bc%vals,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)

          ALLOCATE( bc%vals(nData,ijBeg:ijEnd),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          DO n = 1, nData
            bc%vals(n,:) = holdvals(n)
          ENDDO

          DEALLOCATE ( holdvals,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)

        ENDIF

        tbc%mean = bc%vals(var,:)

      ENDIF    ! active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

  IF (ALLOCATED(morevals)) THEN
    DEALLOCATE ( morevals,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
  ENDIF

  GOTO 999

! error handling --------------------------------------------------------------

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__ )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadTbcSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadTbcSection.F90,v $
! Revision 1.6  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/08/07 16:32:08  rfiedler
! Set bcCode to BC_INFLOW for all inflow BCs so bcCode+BC_RANGE has desired value.
!
! Revision 1.3  2006/08/19 15:38:24  mparmar
! Renamed patch variables
!
! Revision 1.2  2005/04/28 05:45:56  wasistho
! added velocity based inflow BC
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.3  2003/06/04 20:06:38  jferry
! added capability for PEUL BCs with TBCs
!
! Revision 1.2  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.7  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.6  2002/09/28 14:33:24  jferry
! added relative velocity to injection BC
!
! Revision 1.5  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.4  2002/09/25 18:29:57  jferry
! simplified TBC parameter lists
!
! Revision 1.3  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/18 21:50:49  jferry
! Streamlined inelegant coding
!
! Revision 1.1  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
!******************************************************************************







