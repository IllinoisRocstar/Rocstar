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
! Purpose: Wrapper for setting continuous-phase variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: 
!   1. Split from RFLU_SetVarsWrapper because in RungeKuttaMP need to update
!      continuous-phase variables separately and discrete-phase variables can 
!      only be set once their communication is completed.
!
! ******************************************************************************
!
! $Id: RFLU_SetVarsContWrapper.F90,v 1.4 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVarsContWrapper(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, & 
                               RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: RFLU_SetVars

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_UpdateDependentVars
#endif
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVarsContWrapper.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVarsContWrapper',&
  'RFLU_SetVarsContWrapper.F90')

! ******************************************************************************
! Set variables
! ******************************************************************************

! ==============================================================================
! Mixture. NOTE need species conserved state vector in primitive form
! ==============================================================================

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed
#endif

  CALL RFLU_SetVars(pRegion,icgBeg,icgEnd)

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed
#endif

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_UpdateDependentVars(pRegion,icgBeg,icgEnd)
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVarsContWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVarsContWrapper.F90,v $
! Revision 1.4  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2005/11/11 13:26:43  fnajjar
! Initial revision
!
! ******************************************************************************







