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
! Purpose: Initialize user input parameters for mixture and base solver to
!   default values.
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes:
!   1. Variables stored in global type are initialized in RFLU_InitGlobal.
!
! ******************************************************************************
!
! $Id: RFLU_InitInputValues.F90,v 1.36 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitInputValues(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iReg
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitInputValues.F90,v $ $Revision: 1.36 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_InitInputValues',&
  'RFLU_InitInputValues.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing input variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Region-independent (global) values are not initialized here, but in
! RFLU_InitGlobal.F90.
! ******************************************************************************


! ******************************************************************************
! Region-dependent values
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pMixtInput => regions(iReg)%mixtInput

! ==============================================================================
!   Frozen flag
! ==============================================================================

    pMixtInput%frozenFlag = .FALSE.

! ==============================================================================
!   Dimensionality
! ==============================================================================

    pMixtInput%dimens = 3

! ==============================================================================
!   Fluid model
! ==============================================================================

    pMixtInput%fluidModel = FLUID_MODEL_COMP

! ==============================================================================
!   Flow and gas model
! ==============================================================================

    pMixtInput%flowModel = FLOW_EULER
    pMixtInput%gasModel  = GAS_MODEL_TCPERF
    pMixtInput%computeTv = .FALSE.

! ==============================================================================
!   Multi-physics modules
! ==============================================================================

    pMixtInput%turbModel = TURB_MODEL_NONE
    pMixtInput%radiUsed  = .FALSE.  ! No radiation

    pMixtInput%indMfMixt = 0
    pMixtInput%indSd     = 0

! ==============================================================================
!   Spatial discretization
! ==============================================================================

    pMixtInput%spaceDiscr = DISCR_UPW_ROE
    pMixtInput%spaceOrder = DISCR_ORDER_1
    
    pMixtInput%reconst = RECONST_WENO_XYZ

    pMixtInput%cReconstCells = CONSTR_NONE
    pMixtInput%cReconstFaces = CONSTR_NONE

    pMixtInput%cReconstCellsWeight = 1.0_RFREAL
    pMixtInput%cReconstFacesWeight = 1.0_RFREAL

    pMixtInput%stencilDimensCells  = 3
    pMixtInput%stencilDimensFaces  = 3
    pMixtInput%stencilDimensBFaces = 3    
    pMixtInput%spaceOrderBFaces    = DISCR_ORDER_1    

    pMixtInput%timeScheme = TST_HYB5RK
    pMixtInput%cfl        = 3.0_RFREAL
    pMixtInput%epsentr    = 0.05_RFREAL
    pMixtInput%dissFact   = 1.0_RFREAL

! ==============================================================================
!   In-cell test tolerance
! ==============================================================================
    
    pMixtInput%tolerICT = 1.0E-11_RFREAL  
    
! ==============================================================================
!   Grid motion
! ==============================================================================

#ifndef GENX
    pMixtInput%moveGrid      = .FALSE.
    pMixtInput%moveGridType  = MOVEGRID_TYPE_DISP
    pMixtInput%moveGridNIter = 4
    pMixtInput%moveGridSFact = 0.25_RFREAL
#else
    pMixtInput%moveGrid      = .TRUE.
    pMixtInput%moveGridType  = MOVEGRID_TYPE_GENX   
#endif

! ==============================================================================
!   Viscosity model. NOTE reference viscosity set to crazy value so can detect
!   whether it was set in RFLU_DerivedInputValues. If it was not set by user, it
!   is set equal to reference viscosity from REFERENCE input section.
! ==============================================================================

    pMixtInput%viscModel = VISC_SUTHR
    pMixtInput%refTemp   = 110.0_RFREAL
    pMixtInput%refVisc   = REAL(CRAZY_VALUE_INT,RFREAL)
    pMixtInput%suthCoef  = 288.16_RFREAL
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing input variables done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitInputValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitInputValues.F90,v $
! Revision 1.36  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.35  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.34  2007/03/19 21:39:40  haselbac
! Removed init of pMixtInput%nPv
!
! Revision 1.33  2006/10/20 21:24:48  mparmar
! Initialized spaceOrderBFaces with parameter
!
! Revision 1.32  2006/08/19 15:38:44  mparmar
! Added initialization of mixtInput%spaceOrderBFaces
!
! Revision 1.31  2006/04/07 14:43:20  haselbac
! Added init of stencilDimens params
!
! Revision 1.30  2006/01/06 22:06:44  haselbac
! Added seeting of stencilDimens
!
! Revision 1.29  2005/12/25 15:23:16  haselbac
! Added init for constrained reconstruction
!
! Revision 1.28  2005/12/24 21:26:03  haselbac
! Added init of ICT tolerance
!
! Revision 1.27  2005/11/10 22:21:55  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.26  2005/10/31 21:09:34  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.25  2005/10/31 19:25:21  haselbac
! Added init of gasModel
!
! Revision 1.24  2005/10/27 19:43:16  haselbac
! Change init of constr variable
!
! Revision 1.23  2005/10/27 18:56:13  haselbac
! Added init of constr variable
!
! Revision 1.22  2005/08/18 18:47:54  haselbac
! Added init for nPv
!
! Revision 1.21  2005/07/11 19:23:37  mparmar
! Added init of reconst option
!
! Revision 1.20  2005/04/20 14:39:16  haselbac
! Removed CHECK_UNIFLOW code section
!
! Revision 1.19  2005/03/09 14:53:24  haselbac
! Added init of dimensionality
!
! Revision 1.18  2004/11/02 02:28:24  haselbac
! Added initialization for fluid model
!
! Revision 1.17  2004/10/19 19:24:25  haselbac
! Modified initialization of grid-motion input parameters
!
! Revision 1.16  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.15  2004/07/28 15:39:51  jferry
! created global variable for spec use
!
! Revision 1.14  2004/07/08 02:16:52  haselbac
! Added initialization for dissFact
!
! Revision 1.13  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.12  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.11  2004/01/29 22:56:34  haselbac
! Added initialization of indMfMixt
!
! Revision 1.10  2003/12/04 03:23:55  haselbac
! Moved init of global params, changed init of viscosity params
!
! Revision 1.9  2003/11/25 21:02:48  haselbac
! Added initialization of specUsed
!
! Revision 1.8  2003/11/21 22:35:50  fnajjar
! Update Random Number Generator
!
! Revision 1.7  2003/08/20 20:40:07  haselbac
! Bug fix: END DO in wrong place...
!
! Revision 1.6  2003/07/22 01:56:18  haselbac
! Cosmetics only
!
! Revision 1.5  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.4  2003/06/04 22:02:39  haselbac
! Added setting of uniform state (for rfluprep)
!
! Revision 1.3  2003/04/11 20:05:02  fnajjar
! Added default values for viscosity model based on Sutherland Law
!
! Revision 1.2  2003/03/31 16:12:02  haselbac
! Cosmetics, added grid-motion type init
!
! Revision 1.1  2003/01/28 15:53:32  haselbac
! Initial revision, moved from rocflu, use LBOUND and UBOUND
!
! Revision 1.9  2002/11/02 02:04:25  wasistho
! Added TURB statistics
!
! Revision 1.8  2002/10/27 19:13:33  haselbac
! Added setting of moveGrid
!
! Revision 1.7  2002/10/05 19:23:19  haselbac
! GENX integration
!
! Revision 1.6  2002/09/09 15:51:56  haselbac
! global and mixtInput now under region
!
! Revision 1.5  2002/07/25 14:28:10  haselbac
! Added MASTERPROC distinction for output
!
! Revision 1.4  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.3  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.2  2002/05/04 17:10:44  haselbac
! Cosmetic changes
!
! Revision 1.1  2002/03/26 19:25:32  haselbac
! Initial revision
!
! ******************************************************************************







