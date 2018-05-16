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
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DeallocateMemoryWrapper.F90,v 1.16 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DeallocateMemoryWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif

  USE ModInterfaces, ONLY: RFLU_DeallocateMemory, &
                           RFLU_DestroyGrid

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_DeallocMemSol, &
                                     PLAG_RFLU_DeallocMemSolTile, &
                                     PLAG_RFLU_DeallocMemTStep, &
                                     PLAG_RFLU_DeallocMemTStepTile,&
                                     PLAG_INRT_DeallocMemTStep
#endif

#ifdef PERI
  USE ModInterfacesPeriodic, ONLY: PERI_RFLU_DeallocateMemory
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_DeallocateMemory, & 
                                  SPEC_RFLU_DeallocateMemoryEEv
#endif

#ifdef TURB
  USE ModInterfacesTurbulence, ONLY: TURB_RFLU_DeallocateMemory
#endif

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
#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = &
    '$RCSfile: RFLU_DeallocateMemoryWrapper.F90,v $ $Revision: 1.16 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryWrapper',&
  'RFLU_DeallocateMemoryWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory 
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_DeallocateMemory(pRegion)

! ******************************************************************************
! Physical modules
! ******************************************************************************

#ifdef PLAG
! ==============================================================================
! Particles
! ==============================================================================

  IF ( global%plagUsed  .EQV. .TRUE. ) THEN
    pPlag => pRegion%plag
    CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)
    CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
    CALL PLAG_RFLU_DeallocMemTStep(pRegion,pPlag)
    CALL PLAG_RFLU_DeallocMemTStepTile(pRegion)
    CALL PLAG_INRT_DeallocMemTStep(pRegion,pPlag)
  END IF ! plagUsed
#endif

#ifdef RADI
! ==============================================================================
! Radiation
! ==============================================================================

  CALL RADI_DeallocateMemory(pRegion)
#endif

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_DeallocateMemory(pRegion)
    CALL SPEC_RFLU_DeallocateMemoryEEv(pRegion)
  END IF ! global%specUsed
#endif

#ifdef TURB
! ==============================================================================
! Turbulence
! ==============================================================================

  IF ( (pRegion%mixtInput%flowModel == FLOW_NAVST) .AND. &
       (pRegion%mixtInput%turbModel /= TURB_MODEL_NONE) ) THEN
    CALL TURB_RFLU_DeallocateMemory(pRegion)
  END IF ! pRegion%mixtInput%flowModel
#endif

#ifdef PERI
  IF ( pRegion%periInput%flowKind /= OFF ) THEN
    CALL PERI_RFLU_DeallocateMemory(pRegion)
  END IF ! pRegion%periInput%flowKind
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DeallocateMemoryWrapper.F90,v $
! Revision 1.16  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2005/11/27 01:52:19  haselbac
! Added deallocate call for EEv
!
! Revision 1.13  2004/10/19 19:29:11  haselbac
! Adapted to GENX changes, no longer create grid
!
! Revision 1.12  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.11  2004/07/26 19:02:24  fnajjar
! Included call to PLAG_INRT_DeallocMemTStep routine
!
! Revision 1.10  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.9  2004/06/17 23:06:20  wasistho
! added memory deallocation for rocperi
!
! Revision 1.8  2004/03/27 03:03:24  wasistho
! added TURB_RFLU_DeallocateMemory
!
! Revision 1.7  2004/03/19 21:21:30  haselbac
! Cosmetics only
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/02/26 21:02:08  haselbac
! Added PLAG support
!
! Revision 1.4  2004/02/02 22:51:25  haselbac
! Commented out PLAG_DeallocateMemory - temporary measure
!
! Revision 1.3  2003/11/25 21:04:39  haselbac
! Added call to SPEC_RFLU_DeallocateMemory
!
! Revision 1.2  2003/03/15 18:35:03  haselbac
! Adaptation for parallel computations
!
! Revision 1.1  2003/01/28 14:59:34  haselbac
! Initial revision
!
! ******************************************************************************







