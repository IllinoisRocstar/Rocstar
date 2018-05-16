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
! Purpose: Wrapper for setting discrete-phase variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
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
! $Id: RFLU_SetVarsDiscWrapper.F90,v 1.3 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVarsDiscWrapper(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, & 
                               RFLU_ConvertCvPrim2Cons

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_NonCvUpdate
#endif
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVarsDiscWrapper.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVarsDiscWrapper',&
  'RFLU_SetVarsDiscWrapper.F90')

! ******************************************************************************
! Set variables
! ******************************************************************************

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
    CALL PLAG_NonCvUpdate(pRegion)
    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
  END IF ! plagUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVarsDiscWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVarsDiscWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/11/11 13:26:43  fnajjar
! Initial revision
!
! ******************************************************************************







