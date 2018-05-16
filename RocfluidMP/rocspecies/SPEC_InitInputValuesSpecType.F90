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
! Purpose: Initialize user input for species to default values.
!
! Description: None.
!
! Input:
!   regions        Region data
!   iSpecType        Species type
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: SPEC_InitInputValuesSpecType.F90,v 1.8 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_InitInputValuesSpecType(region,iSpecType)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModSpecies, ONLY: t_spec_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iSpecType
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = &
    '$RCSfile: SPEC_InitInputValuesSpecType.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction(global,'SPEC_InitInputValuesSpecType',&
  'SPEC_InitInputValuesSpecType.F90')

! ******************************************************************************
! Initialize
! ******************************************************************************

! ==============================================================================
! Common to all species
! ==============================================================================

  region%specInput%specType(iSpecType)%frozenFlag     = .FALSE.
  region%specInput%specType(iSpecType)%initVal        = 1.0_RFREAL
  region%specInput%specType(iSpecType)%sourceType     = SPEC_SOURCE_TYPE_NONE
  region%specInput%specType(iSpecType)%schmidtNumber  = 1.0_RFREAL

! ==============================================================================
! Specific to species representing particles
! ==============================================================================
  
  region%specInput%specType(iSpecType)%iSpec2iSpecEEv = CRAZY_VALUE_INT
  region%specInput%specType(iSpecType)%iSpecEEv2iSpec = CRAZY_VALUE_INT  
  region%specInput%specType(iSpecType)%diameter       = 0.0_RFREAL
  region%specInput%specType(iSpecType)%puffFactor     = 1.0_RFREAL
  region%specInput%specType(iSpecType)%velocityMethod = SPEC_METHV_FLUIDVEL
  region%specInput%specType(iSpecType)%settlingFlag   = .FALSE.

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_InitInputValuesSpecType

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_InitInputValuesSpecType.F90,v $
! Revision 1.8  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2005/11/27 01:54:46  haselbac
! Added init for EEv variables, cosmetics
!
! Revision 1.5  2005/11/10 02:35:13  haselbac
! Added init for settlingFlag
!
! Revision 1.4  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.3  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.2  2004/01/29 22:59:34  haselbac
! Added Schmidt number
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
!******************************************************************************







