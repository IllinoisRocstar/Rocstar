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
! Purpose: Memory allocation wrapper for rfluconv.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.

!******************************************************************************
!
! $Id: RFLU_AllocMemSolWrapper.F90,v 1.6 2008/12/06 08:44:55 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_AllocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif

  USE RFLU_ModAllocateMemory, ONLY: RFLU_AllocateMemorySolCv

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_AllocMemSol, &
                                     PLAG_RFLU_AllocMemSolTile
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_AllocateMemorySol
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = &
    '$RCSfile: RFLU_AllocMemSolWrapper.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocMemSolWrapper', &
                        'RFLU_AllocMemSolWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  CALL RFLU_AllocateMemorySolCv(pRegion)

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pPlag => pRegion%plag
    CALL PLAG_RFLU_AllocMemSol(pRegion,pPlag)
    CALL PLAG_RFLU_AllocMemSolTile(pRegion)
  END IF ! plagUsed
#endif

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_AllocateMemorySol(pRegion)
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocMemSolWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocMemSolWrapper.F90,v $
! Revision 1.6  2008/12/06 08:44:55  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.3  2004/03/19 21:23:20  haselbac
! Changed allocation logic
!
! Revision 1.2  2004/03/05 22:09:05  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/02/26 21:00:59  haselbac
! Initial revision
!
!******************************************************************************







