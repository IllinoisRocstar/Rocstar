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
! Purpose: Deallocate memory wrapper for rflupost.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.

! ******************************************************************************
!
! $Id: RFLU_DeallocMemSolWrapper.F90,v 1.10 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DeallocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif

  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemorySolCv, & 
                                      RFLU_DeallocateMemorySolDv, &
                                      RFLU_DeallocateMemorySolGv, & 
                                      RFLU_DeallocateMemorySolTv

#ifdef PLAG
  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemorySolTv
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_DeallocMemSol, &
                                     PLAG_RFLU_DeallocMemSolTile
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_DeallocateMemoryEEv, & 
                                  SPEC_RFLU_DeallocateMemorySol
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
    '$RCSfile: RFLU_DeallocMemSolWrapper.F90,v $ $Revision: 1.10 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocMemSolWrapper', &
                        'RFLU_DeallocMemSolWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  CALL RFLU_DeallocateMemorySolCv(pRegion)
  CALL RFLU_DeallocateMemorySolDv(pRegion)
  CALL RFLU_DeallocateMemorySolGv(pRegion)  

  IF ( pRegion%mixtInput%computeTv .EQV. .TRUE. ) THEN 
    CALL RFLU_DeallocateMemorySolTv(pRegion)
  END IF ! pRegion%mixtInput%computeTv

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pPlag => pRegion%plag

    CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)
    CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
  END IF ! plagUsed
#endif

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_DeallocateMemorySol(pRegion)

    IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
      CALL SPEC_RFLU_DeallocateMemoryEEv(pRegion)
    END IF ! pRegion%specInput%nSpeciesEE
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocMemSolWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DeallocMemSolWrapper.F90,v $
! Revision 1.10  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2005/11/27 01:58:43  haselbac
! Added deallocation of EEv
!
! Revision 1.7  2005/11/17 14:49:39  haselbac
! Bug fix: Tv needed by Rocspecies with EE
!
! Revision 1.6  2005/11/11 23:45:59  fnajjar
! Bug fix: Added dealloc of tv for PLAG
!
! Revision 1.5  2005/11/10 02:45:54  haselbac
! Change logic bcos of variable properties
!
! Revision 1.4  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.3  2004/03/19 21:32:14  haselbac
! Changed memory deallocation logic
!
! Revision 1.2  2004/03/05 22:09:05  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/02/26 21:01:22  haselbac
! Initial revision
!
! ******************************************************************************







