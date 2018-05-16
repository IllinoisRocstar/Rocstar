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
! Purpose: Allocate memory wrapper.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: 
!   1. Do not allocate memory for particles because treated separately.
!
! ******************************************************************************
!
! $Id: RFLU_AllocMemSolWrapper.F90,v 1.6 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_AllocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModAllocateMemory, ONLY: RFLU_AllocateMemorySolCv, & 
                                    RFLU_AllocateMemorySolDv, & 
                                    RFLU_AllocateMemorySolGv

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_AllocateMemoryEEv, &
                                  SPEC_RFLU_AllocateMemorySol
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

  RCSIdentString = '$RCSfile: RFLU_AllocMemSolWrapper.F90,v $ $Revision: 1.6 $'

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

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_AllocateMemorySolCv(pRegion)
  CALL RFLU_AllocateMemorySolGv(pRegion)  
  CALL RFLU_AllocateMemorySolDv(pRegion)

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_AllocateMemorySol(pRegion)
    
    IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
      CALL SPEC_RFLU_AllocateMemoryEEv(pRegion)
    END IF ! pRegion%specInput%nSpeciesEE    
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocMemSolWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocMemSolWrapper.F90,v $
! Revision 1.6  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/24 02:14:27  haselbac
! Bug fix: Allocate dv bcos of NSCBC code
!
! Revision 1.3  2005/11/27 01:57:47  haselbac
! Added allocation of EEv
!
! Revision 1.2  2005/11/10 02:41:02  haselbac
! Now alloc gv outside GENX bcos of variable properties
!
! Revision 1.1  2005/04/15 15:08:11  haselbac
! Initial revision
!
! ******************************************************************************







