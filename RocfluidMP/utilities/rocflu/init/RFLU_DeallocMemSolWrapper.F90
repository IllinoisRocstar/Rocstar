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
! Purpose: Deallocate memory wrapper.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: 
!   1. Do not deallocate memory for particles because treated separately.
!
! ******************************************************************************
!
! $Id: RFLU_DeallocMemSolWrapper.F90,v 1.6 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DeallocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemorySolCv, & 
                                      RFLU_DeallocateMemorySolDv, & 
                                      RFLU_DeallocateMemorySolGv

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_DeallocateMemoryEEv, & 
                                  SPEC_RFLU_DeallocateMemorySol
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

  RCSIdentString = '$RCSfile: RFLU_DeallocMemSolWrapper.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocMemSolWrapper', &
                        'RFLU_DeallocMemSolWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_DeallocateMemorySolCv(pRegion)
  CALL RFLU_DeallocateMemorySolGv(pRegion)  
  CALL RFLU_DeallocateMemorySolDv(pRegion)  

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_DeallocateMemorySol(pRegion)
    
    IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
      CALL SPEC_RFLU_DeallocateMemoryEEv(pRegion)
    END IF ! pRegion%specInput%nSpeciesEE      
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

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
! Revision 1.6  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/24 02:14:48  haselbac
! Bug fix: Deallocate dv bcos of NSCBC code
!
! Revision 1.3  2005/11/27 01:58:02  haselbac
! Added deallocation of EEv
!
! Revision 1.2  2005/11/10 02:41:37  haselbac
! Now dealloc gv outside GENX bcos of variable properties
!
! Revision 1.1  2005/04/15 15:08:12  haselbac
! Initial revision
!
! ******************************************************************************







