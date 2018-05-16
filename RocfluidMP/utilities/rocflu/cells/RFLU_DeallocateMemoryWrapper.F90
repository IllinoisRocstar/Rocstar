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
! Purpose: Deallocate memory wrapper.
!
! Description: None.
!
! Input:
!   pRegion	Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_DeallocateMemoryWrapper.F90,v 1.3 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_DeallocateMemoryWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI
  
  USE RFLU_ModCellMapping
  USE RFLU_ModFaceList
    
  USE ModInterfaces, ONLY: RFLU_DeallocateMemory,RFLU_DestroyGrid
  
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

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_DeallocateMemoryWrapper.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryWrapper', &
                        'RFLU_DeallocateMemoryWrapper.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Destroy data structure
! ******************************************************************************

  CALL RFLU_DestroyCellMapping(pRegion)
  CALL RFLU_DestroyFaceList(pRegion,DESTROY_FACE_ALL)

! ******************************************************************************
! Destroy grid
! ******************************************************************************

  CALL RFLU_DestroyGrid(pRegion)

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DeallocateMemoryWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
!******************************************************************************







