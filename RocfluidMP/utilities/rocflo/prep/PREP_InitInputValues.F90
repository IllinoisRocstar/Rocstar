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
! Purpose: initialize user input parameters to default values.
!
! Description: none.
!
! Input: global, regions
!
! Output: global, regions
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_InitInputValues.F90,v 1.5 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE InitInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModError
  USE ModParameters
  USE PREP_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_mixt_input), POINTER :: input
  TYPE(t_global), POINTER     :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'InitInputValues',&
  'PREP_InitInputValues.F90' )

! global values ---------------------------------------------------------------
! level, formats and flow-type

  global%startLevel  = 1 
  global%gridFormat  = FORMAT_ASCII 
  global%solutFormat = FORMAT_ASCII 
  global%flowType    = FLOW_UNSTEADY

! reference values

  global%refVelocity = 100.
  global%refPressure = 1.E+5
  global%refDensity  = 1.2
  global%refCp       = 1004.5
  global%refGamma    = 1.4
  global%refLength   = 1.0
  global%refREnum    = 100.
  global%prLam       = 0.72
  global%prTurb      = 0.9
  global%scnLam      = 0.22
  global%scnTurb     = 0.9

! multiphysics modules

  global%peulUsed = .false.
  global%plagUsed = .false.
  global%inrtUsed = .false.

! region related values -------------------------------------------------------

  DO iReg=1,global%nRegions

    regions(iReg)%procid    = global%myProcid
    regions(iReg)%nDumCells = 2

! - input parameters

    input => regions(iReg)%mixtInput

    input%flowModel   = FLOW_EULER
    input%turbModel   = TURB_MODEL_NONE
    input%moveGrid    = .false.
    input%computeTv   = .false.

    input%prepIniCase = INITFLO_UNIFORM
    input%gasModel   = GAS_MODEL_TCPERF

    input%radiUsed    = .false.  ! no radiation

    input%iniVelX     = 100._RFREAL
    input%iniVelY     = 0._RFREAL
    input%iniVelZ     = 0._RFREAL
    input%iniPress    = 1.E+5_RFREAL
    input%iniDens     = 1.2_RFREAL
    input%iniXsplit   = 0._RFREAL

  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_InitInputValues.F90,v $
! Revision 1.5  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/10/31 21:09:38  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.2  2005/09/21 20:35:54  wasistho
! initialize iniCase and iniXsplit
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.8  2004/07/27 03:34:30  wasistho
! add additional default input variables
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/03 23:55:41  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.5  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/20 19:35:43  haselbac
! Modified RegFun call to avoid probs with long 'PREP_InitInputValues.F90' names
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
!******************************************************************************








