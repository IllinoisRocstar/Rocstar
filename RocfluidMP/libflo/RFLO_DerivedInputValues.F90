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
! Purpose: set values derived from user input.
!
! Description: none.
!
! Input: regions = input parameters for all regions.
!
! Output: regions = numerical parameters and no. of equations for each region.
!
! Notes: dimensions of work arrays are also set here (for the mixture).
!
!******************************************************************************
!
! $Id: RFLO_DerivedInputValues.F90,v 1.8 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_DerivedInputValues( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModInterfaces, ONLY : RFLO_SetMstageCoeffs
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev

! ... local variables
  INTEGER :: dim1d
  LOGICAL :: someMoved

  TYPE(t_mixt_input), POINTER :: input
  TYPE(t_global), POINTER     :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_DerivedInputValues',&
  'RFLO_DerivedInputValues.F90' )

! global values ---------------------------------------------------------------

  IF (global%refREnum > 0._RFREAL) THEN
    global%refVisc = global%refDensity*global%RefVelocity*global%refLength/ &
                     global%refREnum
  ELSE
    global%refVisc = 0._RFREAL
  ENDIF

! compute cell and face centers?

  someMoved = .false.
  DO iReg=1,global%nRegions
    IF (regions(iReg)%mixtInput%moveGrid .eqv. .true.) THEN
      someMoved = .true.
      EXIT
    ENDIF
  ENDDO

  IF ( (someMoved.eqv..true.) .AND. &
      (global%moveGridScheme==MOVEGRID_VMS .OR. &
       global%moveGridScheme==MOVEGRID_FOMS)) THEN
    global%calcCellCtr = .TRUE.
    global%calcFaceCtr = .TRUE.
  ENDIF

  IF (global%aeroCoeffs==ACTIVE) THEN
    global%calcFaceCtr = .TRUE.
  ENDIF

! region related data ---------------------------------------------------------

  DO iReg=1,global%nRegions

    input => regions(iReg)%mixtInput

! - number of dummy cells

    IF (input%spaceOrder > DISCR_ORDER_2) THEN     ! > 2nd-order
      regions(iReg)%nDumCells = 3                  ! (LES may add more)
    ENDIF

! - determine if we need to compute transport variables

    input%computeTv = (input%flowModel == FLOW_NAVST)

! - mixture constants, no. of variables, indexing of Cp and Mol vectors

    DO iLev=1,regions(iReg)%nGridLevels
      regions(iReg)%levels(iLev)%mixt%nDv     = 6   ! dependent variables
      IF (input%computeTv.eqv..true.) THEN
        regions(iReg)%levels(iLev)%mixt%nTv   = 2   ! transport variables
      ELSE
        regions(iReg)%levels(iLev)%mixt%nTv   = 0   ! transport variables
      ENDIF
      regions(iReg)%levels(iLev)%mixt%nGv     = 2   ! gas variables
      IF (input%flowModel == FLOW_NAVST) THEN
        regions(iReg)%levels(iLev)%mixt%nGrad = 12  ! gradients
      ELSE
        regions(iReg)%levels(iLev)%mixt%nGrad = 0   ! gradients
      ENDIF
      regions(iReg)%levels(iLev)%mixt%indCp   = 0   ! pointer to Cp
      regions(iReg)%levels(iLev)%mixt%indMol  = 0   ! pointer to Mol

      regions(iReg)%levels(iLev)%mixt%prLam   = global%prLam
      regions(iReg)%levels(iLev)%mixt%prTurb  = global%prTurb
      regions(iReg)%levels(iLev)%mixt%scnLam  = global%scnLam
      regions(iReg)%levels(iLev)%mixt%scnTurb = global%scnTurb
    ENDDO

    IF (input%refVisc < 0._RFREAL) input%refVisc = global%refVisc

! - grid motion

    IF (.NOT. (input%moveGrid.eqv..true.)) THEN
      global%moveGridScheme = MOVEGRID_BLOCKS
    ENDIF

! - grid variables (grid speeds)

    DO iLev=1,regions(iReg)%nGridLevels
      IF (input%moveGrid.eqv..true.) THEN
        regions(iReg)%levels(iLev)%grid%indSvel = 1
      ELSE
        regions(iReg)%levels(iLev)%grid%indSvel = 0
      ENDIF
    ENDDO

! - time-stepping scheme

    IF (global%flowType == FLOW_STEADY) THEN     ! steady flow
      input%timeScheme  = TST_HYB5RK
      global%tstepOrder = 0
    ELSE                                         ! unsteady flow
      IF (global%solverType == SOLV_EXPLICIT) THEN
        input%timeScheme  = TST_STD4RK
        global%tstepOrder = 4
      ELSE
        input%timeScheme  = TST_HYB5RK
        global%startLevel = 1     ! no FMG here
      ENDIF
    ENDIF

    CALL RFLO_SetMstageCoeffs( global,regions(iReg)%mixtInput,global%nrkSteps )

! - dimensions of work arrays

    dim1d = MAX(regions(iReg)%levels(1)%grid%ipc, &
                regions(iReg)%levels(1)%grid%jpc, &
                regions(iReg)%levels(1)%grid%kpc) + 6  ! 2*3 dummy cells
    regions(iReg)%dimWork1D    = dim1D*dim1D
    regions(iReg)%dimWork2D(:) = 1

  ENDDO   ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_DerivedInputValues.F90,v $
! Revision 1.8  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/10/23 18:20:53  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.5  2006/03/24 05:01:27  wasistho
! calc faceCtr for active global%aeroCoeffs
!
! Revision 1.4  2005/10/20 06:52:16  wasistho
! assign value to calcFaceCtr
!
! Revision 1.3  2005/10/17 22:33:40  wasistho
! moved calcCellCtr test to from InitFlowSolver to DerivedInputValues
!
! Revision 1.2  2005/06/26 08:42:50  wasistho
! set moveGridScheme to type 0 if no grid motion
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.30  2004/11/17 16:12:34  haselbac
! Added global pointer to RFLO_SetMStageCoeffs interface
!
! Revision 1.29  2004/04/08 03:15:57  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.28  2004/03/05 22:08:59  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.27  2004/03/03 23:55:38  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.26  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.23  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.22  2003/05/20 22:16:16  jblazek
! Corrected bug in viscosity model input.
!
! Revision 1.21  2003/05/20 21:13:20  jblazek
! Reworked RK coefficients.
!
! Revision 1.20  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.19  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.18  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.17  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.16  2002/08/30 19:08:58  jblazek
! Dimensions of work arrays now set in derivedInputValues.
!
! Revision 1.15  2002/08/29 23:37:52  jblazek
! Changed name of index pointer for grid speeds. Set nTv=2 for N.-S. equs.
!
! Revision 1.14  2002/08/21 22:47:49  wasistho
! Set nTv = 4 if turbulence is active
!
! Revision 1.13  2002/06/30 00:02:50  jblazek
! Set nTv=2 for viscous flows.
!
! Revision 1.12  2002/05/21 01:50:16  wasistho
! add viscous terms
!
! Revision 1.11  2002/03/30 00:50:48  jblazek
! Cleaned up with flint.
!
! Revision 1.10  2002/03/18 21:56:39  jblazek
! Finished multiblock and MPI.
!
! Revision 1.9  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.8  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
! Revision 1.7  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/06 00:15:39  jblazek
! Improved injection BC. Added pointers to gradients.
!
! Revision 1.5  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.4  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.3  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/12 00:02:48  jblazek
! Added postprocessor.
!
!******************************************************************************







