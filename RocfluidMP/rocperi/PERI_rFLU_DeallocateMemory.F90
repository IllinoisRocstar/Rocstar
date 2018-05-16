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
! Purpose: Deallocate memory for rocperi.
!
! Description: None.
!
! Input: region = Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PERI_rFLU_DeallocateMemory.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PERI_RFLU_DeallocateMemory(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

!  USE PERI_ModInterfaces, ONLY: PERI_FluDeallocateMemorySol, & 
!                                PERI_FluDeallocateMemoryTStep

  IMPLICIT NONE

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: region

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PERI_rFLU_DeallocateMemory.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction(global,'PERI_RFLU_DeallocateMemory',&
  'PERI_rFLU_DeallocateMemory.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory for rocperi...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

!  CALL PERI_FluDeallocateMemorySol(region)
!  CALL PERI_FluDeallocateMemoryTStep(region)

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                 'Deallocating memory for rocperi done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PERI_RFLU_DeallocateMemory

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_rFLU_DeallocateMemory.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/06/17 23:09:52  wasistho
! initial import PERI_RFLU_DeallocateMemory
!
!
!
! ******************************************************************************







