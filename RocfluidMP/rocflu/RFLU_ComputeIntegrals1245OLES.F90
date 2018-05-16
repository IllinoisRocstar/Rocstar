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
! Purpose: Compute integrals 1,2,4,5 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to current region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ComputeIntegrals1245OLES.F90,v 1.5 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegrals1245OLES(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_ComputeIntegral1OLES, & 
                           RFLU_ComputeIntegral2OLES, & 
                           RFLU_ComputeIntegral4OLES, & 
                           RFLU_ComputeIntegral5OLES

  IMPLICIT NONE

! --- parameters      

  TYPE(t_region), POINTER :: pRegion

! --- local variables

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Start
! ==============================================================================

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegrals1245OLES',&
  'RFLU_ComputeIntegrals1245OLES.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                             'integrals...'
  END IF ! global%verbLevel

! ==============================================================================
! Compute integrals
! ==============================================================================

  CALL RFLU_ComputeIntegral1OLES(pRegion)
  CALL RFLU_ComputeIntegral2OLES(pRegion)
  CALL RFLU_ComputeIntegral4OLES(pRegion)
  CALL RFLU_ComputeIntegral5OLES(pRegion)

! ==============================================================================
! End
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                             'integrals done.'
  END IF ! global%verbLevel  

  CALL DeregisterFunction(global)     

END SUBROUTINE RFLU_ComputeIntegrals1245OLES







