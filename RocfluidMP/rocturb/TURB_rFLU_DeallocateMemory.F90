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
! Purpose: Deallocate memory for turbulence.
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
! $Id: TURB_rFLU_DeallocateMemory.F90,v 1.4 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_DeallocateMemory(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

!  USE TURB_ModInterfaces, ONLY: TURB_FluDeallocateMemorySol, & 
!                                TURB_FluDeallocateMemoryTStep

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

  RCSIdentString = '$RCSfile: TURB_rFLU_DeallocateMemory.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction(global,'TURB_RFLU_DeallocateMemory',&
  'TURB_rFLU_DeallocateMemory.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory for turbulence...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

!  CALL TURB_FluDeallocateMemorySol(region)
!  CALL TURB_FluDeallocateMemoryTStep(region)

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                 'Deallocating memory for turbulence done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_RFLU_DeallocateMemory

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLU_DeallocateMemory.F90,v $
! Revision 1.4  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/17 23:09:19  wasistho
! change interface module to TURB_ModInterfaces
!
! Revision 1.1  2004/03/27 02:19:15  wasistho
! added routines specific for Rocflu
!
!
! ******************************************************************************







