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
! Purpose: Initialize user input for species to default values.
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_InitInputValues.F90,v 1.7 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_InitInputValues(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModSpecies, ONLY: t_spec_input
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
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_InitInputValues.F90,v $ $Revision: 1.7 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_InitInputValues',&
  'SPEC_InitInputValues.F90')

! ******************************************************************************
! Initialize
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    regions(iReg)%specInput%usedFlag   = .FALSE.
    regions(iReg)%specInput%nSpecies   = 0
    regions(iReg)%specInput%nSpeciesEE = 0    
    regions(iReg)%specInput%sourceFlag = .FALSE.
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_InitInputValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_InitInputValues.F90,v $
! Revision 1.7  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/11/27 01:54:17  haselbac
! Added init for nSpeciesEE, cosmetics
!
! Revision 1.4  2005/04/20 14:44:36  haselbac
! Removed setting of unifSpec
!
! Revision 1.3  2004/07/28 15:31:34  jferry
! added USED field to SPECIES input section
!
! Revision 1.2  2004/04/01 21:31:37  haselbac
! Added setting of sourceFlag
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







