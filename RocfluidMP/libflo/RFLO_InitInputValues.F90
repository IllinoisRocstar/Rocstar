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
! Purpose: initialize user input parameters for mixture & base solver
!          to default values.
!
! Description: none.
!
! Input: none.
!
! Output: regions = initial input values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InitInputValues.F90,v 1.26 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InitInputValues( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  TYPE(t_mixt_input), POINTER :: input
  TYPE(t_patch), POINTER      :: patch
  TYPE(t_global), POINTER     :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_InitInputValues',&
  'RFLO_InitInputValues.F90' )

! global values ---------------------------------------------------------------
! regions per processor

  global%nRegionsProc = 1

! random number generator

  global%randSeedOffset = 0
  global%randSeedType   = RAND_SEED_TYPE_FIXED

! formats for grid & solution

  global%gridFormat  = FORMAT_ASCII
  global%solutFormat = FORMAT_ASCII

! reference values

  global%refVelocity = 100._RFREAL
  global%refPressure = 1.E+5_RFREAL
  global%refDensity  = 1.2_RFREAL
  global%refCp       = 1004.5_RFREAL
  global%refGamma    = 1.4_RFREAL
  global%refLength   = 1.0_RFREAL
  global%refREnum    = 100._RFREAL
  global%prLam       = 0.72_RFREAL
  global%prTurb      = 0.9_RFREAL
  global%scnLam      = 0.22_RFREAL
  global%scnTurb     = 0.9_RFREAL

! acceleration

  global%accelOn = .false.
  global%accelX  = 0._RFREAL
  global%accelX  = 0._RFREAL
  global%accelX  = 0._RFREAL

! grid motion

  global%moveGridScheme    = MOVEGRID_BLOCKS
  global%moveGridNiter     = 10
  global%moveGridViter     = 5
  global%moveGridSiter     = 20
  global%moveGridWeight    = 2._RFREAL
  global%moveGridAmplifX   = 1._RFREAL
  global%moveGridAmplifY   = 1._RFREAL
  global%moveGridAmplifZ   = 1._RFREAL
  global%moveGridPower     = 1._RFREAL
  global%moveGridNbour     = 12
  global%moveGridNsmatch   = 20
  global%moveGridOrthDir   = OFF
  global%moveGridOrthWghtX = 0._RFREAL
  global%moveGridOrthWghtY = 0._RFREAL
  global%moveGridOrthWghtZ = 0._RFREAL
  global%moveGridOrthCell  = 0._RFREAL

! multigrid

  global%startLevel = 1
  global%cycleType  = MGCYCLE_NO
  global%refineIter = 9999999

! position of probe(s)

  global%nProbes        = 0
  global%probeSaveTime  = 0._RFREAL
  global%probeSaveIter  = 1
  global%probeOpenClose = .false.

! thrust

  global%thrustType      = THRUST_NONE
  global%thrustPlane     = XCOORD
  global%thrustCoord     = 0._RFREAL
  global%thrustPamb      = 1.E+5_RFREAL
  global%thrustSaveTime  = 0._RFREAL
  global%thrustSaveIter  = 1
  global%thrustOpenClose = .false.

! global forces and aerodynamic coefficients

  global%forcesOn = FORCES_NONE
  global%forceX   = 0._RFREAL
  global%forceY   = 0._RFREAL
  global%forceZ   = 0._RFREAL

  global%aeroCoeffs     = OFF
  global%forceRefLength = 1._RFREAL
  global%forceRefArea   = 1._RFREAL
  global%forceRefXCoord = 0._RFREAL
  global%forceRefYCoord = 0._RFREAL
  global%forceRefZCoord = 0._RFREAL
  global%acBndBoxXmin   = -HUGE( 1._RFREAL )
  global%acBndBoxXmax   =  HUGE( 1._RFREAL ) 
  global%acBndBoxYmin   = -HUGE( 1._RFREAL )
  global%acBndBoxYmax   =  HUGE( 1._RFREAL )
  global%acBndBoxZmin   = -HUGE( 1._RFREAL )
  global%acBndBoxZmax   =  HUGE( 1._RFREAL )

! time stepping

  global%flowType      = FLOW_STEADY
  global%solverType    = SOLV_EXPLICIT
  global%currentIter   = 0
  global%maxIter       = 10000
  global%writeIter     = 1000
  global%printIter     = 1
  global%dtImposed     = 1.E-5_RFREAL
#ifndef GENX
  global%timeStamp     = 0.E+0_RFREAL
#endif
  global%maxTime       = 1.E-4_RFREAL
  global%writeTime     = 5.E-5_RFREAL
  global%printTime     = 1.E-15_RFREAL
  global%resTol        = 1.E-5_RFREAL
  global%tstepOrder    = 2
  global%maxSubIter    = 100
  global%tolSubIter    = 1.E-2_RFREAL
  global%predictSol    = .true.
  global%dualTstSource = .false.
  global%dtFixed       = .true.

  global%rkScheme      = RK_SCHEME_4_CLASSICAL

#ifdef STATS
! time averaged statistics

  global%doStat     = 0
  global%reStat     = 0
  global%mixtNStat  = 0
  global%turbNStat  = 0
  global%integrTime = 0._RFREAL
#endif

! boundary conditions

  global%infloNijk    = NIJK_INFLOW_INIT
  global%internDeform = 0

! multiphysics modules

  global%peulUsed = .false.
  global%plagUsed = .false.
  global%inrtUsed = .false.

! post processing

  global%postPlotType  = PLOT_GRID_FLOW 
  global%postOutFmt    = PLOT_FMT_TECASCII
  global%postStatsFlag = .FALSE.
  global%postTurbFlag  = .FALSE.
  global%postPlagFlag  = .FALSE.
  global%postRadiFlag  = .FALSE.
  global%postSpecFlag  = .FALSE.

! region related values -------------------------------------------------------

  DO iReg=1,global%nRegions

    regions(iReg)%nDumCells  = 2
    regions(iReg)%blockShape = REGION_SHAPE_NORMAL

    input => regions(iReg)%mixtInput

    input%flowModel = FLOW_EULER
    input%turbModel = TURB_MODEL_NONE
    input%moveGrid  = .false.
    input%computeTv = .false.

    input%frozenFlag = .FALSE.
    input%gasModel   = GAS_MODEL_TCPERF

    input%radiUsed  = .false.  ! no radiation

    input%spaceDiscr   = DISCR_CEN_SCAL
    input%spaceOrder   = DISCR_ORDER_2
    input%timeScheme   = TST_HYB5RK
    input%cfl          = 3.0_RFREAL
    input%smoocf       = -1.0_RFREAL
    input%vis2         = 0.5_RFREAL
    input%vis4         = 1._RFREAL/128._RFREAL
    input%pSwitchType  = PSWITCH_STD
    input%pSwitchOmega = 0.5_RFREAL
    input%limfac       = 5.0_RFREAL
    input%epsentr      = 0.05_RFREAL
    input%faceEdgeAvg  = FE_AVG_UNIFORM

    input%viscModel    = VISC_SUTHR
    input%refVisc      = global%refDensity*global%refVelocity &
                       * global%refLength/global%refREnum
    input%refTemp      = 110.0_RFREAL
    input%suthCoef     = 288.16_RFREAL

! - patch related values

    DO iPatch=1,regions(iReg)%nPatches
      patch => regions(iReg)%levels(1)%patches(iPatch)   ! only finest level
      patch%globalAeroCoeffs  = .false.
      patch%mixt%setMotion = .false.
      patch%mixt%bndVel(:) = 0._RFREAL
      patch%mixt%amplitude = 0._RFREAL
    ENDDO  ! iPatch

  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InitInputValues.F90,v $
! Revision 1.26  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.25  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.24  2006/08/19 15:38:05  mparmar
! Renamed patch variables
!
! Revision 1.23  2006/05/09 23:34:57  wasistho
! initialized global%dtFixed
!
! Revision 1.22  2006/03/24 05:00:12  wasistho
! huged initial ranges of acBndBox
!
! Revision 1.21  2006/03/18 13:27:17  wasistho
! initialized orthDir and orthWghtX,Y,Z
!
! Revision 1.20  2006/03/14 04:36:32  wasistho
! removed initialization of patch%bndFlat
!
! Revision 1.19  2006/03/10 01:44:44  wasistho
! changed coeffsBox to acBndBox
!
! Revision 1.18  2006/03/10 00:58:07  wasistho
! initialized aerodynamic coeffs data
!
! Revision 1.17  2006/03/08 06:36:15  wasistho
! added siter and viter
!
! Revision 1.16  2006/03/04 04:34:01  wasistho
! initialized region%blockShape
!
! Revision 1.15  2005/12/05 05:56:09  wasistho
! initialize patch%bndFlat
!
! Revision 1.14  2005/12/01 17:06:57  fnajjar
! Added default value to randSeedType
!
! Revision 1.13  2005/12/01 08:58:49  wasistho
! replaced 10000 by NIJK_INFLOW_INIT
!
! Revision 1.12  2005/11/18 07:53:58  wasistho
! initialized valMixt%amplitude
!
! Revision 1.11  2005/11/10 22:19:29  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.10  2005/10/31 21:09:33  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.9  2005/10/28 22:47:53  wasistho
! initialize orthocell
!
! Revision 1.8  2005/09/30 01:02:57  wasistho
! initialize global%internDeform
!
! Revision 1.7  2005/09/20 23:57:20  wasistho
! initialize global%infloNijk
!
! Revision 1.6  2005/09/09 03:23:41  wasistho
! initialize boundary patch motion variables
!
! Revision 1.5  2005/08/28 23:48:32  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.4  2005/08/18 19:48:11  wasistho
! added moveGridNsmatch
!
! Revision 1.3  2005/06/04 01:02:03  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.2  2005/06/02 22:58:48  wasistho
! initiate global%moveGridAmplif and moveGridPower
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.29  2004/11/17 16:12:56  haselbac
! Added initialization for rkScheme
!
! Revision 1.28  2004/09/02 02:34:37  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.27  2004/07/28 01:52:56  wasistho
! initialize global%post... variables
!
! Revision 1.26  2004/04/08 03:15:41  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.25  2004/03/05 22:08:59  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.24  2004/03/03 23:55:38  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.23  2004/03/02 21:49:20  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.22  2003/11/21 22:30:04  fnajjar
! Added global seed offset for Random Number Generator
!
! Revision 1.21  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.20  2003/08/28 20:05:34  jblazek
! Added acceleration terms.
!
! Revision 1.19  2003/08/25 21:51:23  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.18  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.17  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.16  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.15  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.14  2003/04/11 20:02:35  fnajjar
! Added default values for viscosity model based on Sutherland Law
!
! Revision 1.13  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.12  2002/11/02 01:46:53  wasistho
! Added TURB statistics
!
! Revision 1.11  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.10  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.9  2002/07/25 00:39:54  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.8  2002/06/14 20:53:05  wasistho
! add time avg statistics
!
! Revision 1.7  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.6  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
! Revision 1.5  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.3  2002/02/09 01:47:00  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.2  2002/02/01 21:04:25  jblazek
! Streamlined time stepping routine.
!
! Revision 1.1  2002/01/12 00:02:48  jblazek
! Added postprocessor.
!
!******************************************************************************







