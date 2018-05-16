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
! Purpose: Allocate memory wrapper.
!
! Description: None.
!
! Input:
!   pRegion	Region pointer
!   allocMode 	Allocation mode
!
! Output: None.
!
! Notes: 
!   1. At present, do not use allocMode.
!
!******************************************************************************
!
! $Id: RFLU_AllocateMemoryWrapper.F90,v 1.4 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_AllocateMemoryWrapper(pRegion,allocMode)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModParameters
  USE ModMPI
  
  USE RFLU_ModCellMapping
  USE RFLU_ModFaceList
  
  USE ModInterfaces, ONLY: RFLU_AllocateMemory,RFLU_CreateGrid  
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: allocMode
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_AllocateMemoryWrapper.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocateMemoryWrapper', &
                        'RFLU_AllocateMemoryWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory for grid and data structures 
! ******************************************************************************

  CALL RFLU_CreateGrid(pRegion)
  CALL RFLU_CreateCellMapping(pRegion)
  CALL RFLU_CreateFaceList(pRegion)

! ==============================================================================  
! Grid speeds
! ==============================================================================  

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN
    IF ( pGrid%nFaces > 0 ) THEN 
      ALLOCATE(pGrid%gs(pGrid%nFaces),STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%grid%gs')
      END IF ! global%error  

      pGrid%gs(:) = 0.0_RFREAL
    ELSE 
      NULLIFY(pGrid%gs)
    END IF ! pGrid%nFaces   
    
    IF ( pRegion%grid%nPatches > 0 ) THEN
      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)
     
        IF ( pPatch%nBFaces > 0 ) THEN 
          ALLOCATE(pPatch%gs(pPatch%nBFaces),STAT=errorFlag)
          global%error = errorFlag
          IF (global%error /= 0) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%gs')
          END IF ! global%error 

          pPatch%gs(:) = 0.0_RFREAL
        ELSE 
          NULLIFY(pPatch%gs)
        END IF ! pPatch%nBFaces     
      END DO ! iPatch                   
    END IF ! pRegion%grid%nPatches 
  END IF ! pRegion%mixtInput%moveGrid

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocateMemoryWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocateMemoryWrapper.F90,v $
! Revision 1.4  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2003/04/12 22:20:51  haselbac
! Added allocation for grid speeds
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
!******************************************************************************







