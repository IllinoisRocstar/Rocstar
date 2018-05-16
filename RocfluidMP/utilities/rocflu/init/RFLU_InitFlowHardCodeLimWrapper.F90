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
! Purpose: Wrapper for initialize part of flow field in a region using 
!   hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowHardCodeLimWrapper.F90,v 1.3 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowHardCodeLimWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, & 
                               RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: RFLU_InitFlowHardCodeLim, & 
                           RFLU_SetGasVars
                               
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
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
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowHardCodeLimWrapper.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowHardCodeLimWrapper', &
                        'RFLU_InitFlowHardCodeLimWrapper.F90')

#ifdef SPEC
! ******************************************************************************
! Species. NOTE needs to be done first so can define gas properties when these
! are influenced by species. NOTE initialize species state vector in primitive 
! state, so need to convert species state vector, which was already read in, to 
! conservative form.
! ******************************************************************************

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)      
! TEMPORARY  
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
! END TEMPORARY
  END IF ! global%specUsed
#endif

! ******************************************************************************
! Mixture
! ******************************************************************************

  CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)
  CALL RFLU_InitFlowHardCodeLim(pRegion)
  
! ******************************************************************************
! Physical modules
! ******************************************************************************

! TO BE DONE

! END TO BE DONE  
  
#ifdef SPEC
! ******************************************************************************
! Convert species state vector to conservative state
! ******************************************************************************

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed
#endif  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowHardCodeLimWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowHardCodeLimWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/11/10 02:54:32  haselbac
! Initial revision
!
! ******************************************************************************







