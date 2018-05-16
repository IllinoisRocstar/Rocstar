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
! Purpose: Enforce bounds on species.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. At present, only enforce minimum bound. In future, could also enforce
!      maximum bound on each individual component, and that the sum is equal
!      to the mixture density.
!
!******************************************************************************
!
! $Id: SPEC_RFLU_EnforceBounds.F90,v 1.3 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_EnforceBounds(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region 
  USE ModMPI
  
  USE ModInterfaces, ONLY: RFLU_DecidePrint
    
  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,iSpec,negValCntr
  REAL(RFREAL) :: negValMin
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_global), POINTER :: global
  
  RCSIdentString = '$RCSfile: SPEC_RFLU_EnforceBounds.F90,v $ $Revision: 1.3 $'
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_EnforceBounds',&
  'SPEC_RFLU_EnforceBounds.F90')

! ==============================================================================
! Set pointers and variables
! ==============================================================================

  pCvMixt => pRegion%mixt%cv
  pCvSpec => pRegion%spec%cv
  pGrid   => pRegion%grid
  
  negValCntr = 0
  negValMin  = HUGE(1.0_RFREAL)
  
! ******************************************************************************
! Enforce bounds
! ******************************************************************************

  DO icg = 1,pGrid%nCells
    DO iSpec = 1,pRegion%specInput%nSpecies
      IF ( pCvSpec(iSpec,icg) < 0.0_RFREAL ) THEN
        negValCntr = negValCntr + 1 
        negValMin = MIN(negValMin,pCvSpec(iSpec,icg))        
        pCvSpec(iSpec,icg) = 0.0_RFREAL                   
      END IF ! pCvSpec
    END DO ! iSpec
  END DO ! icg
           
! ******************************************************************************
! Write information
! ******************************************************************************
           
  IF ( (global%verbLevel > VERBOSE_LOW) .AND. &
       (global%myProcid == MASTERPROC) ) THEN 
    IF ( (RFLU_DecidePrint(global) .EQV. .TRUE.) .AND. (negValCntr > 0) ) THEN
      global%warnCounter = global%warnCounter + 1
     
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            '*** WARNING *** Detected negative species mass fractions!'
            
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%flowType
      
      WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME,'Runge-Kutta stage:', &
                                     pRegion%irkStep      

      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, & 
            'Number of locations with negative species mass fractions:', &
            negValCntr
      WRITE(STDOUT,'(A,3X,A,1X,E23.16)') SOLVER_NAME, &
            'Largest negative species mass fraction:',negValMin
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
            'Negative species mass fractions set to zero.'            
    END IF ! RFLU_DecidePrint
  END IF ! global%verbLevel           
           
! ******************************************************************************
! End
! ******************************************************************************  
      
  CALL DeregisterFunction(global)  
    
END SUBROUTINE SPEC_RFLU_EnforceBounds


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: SPEC_RFLU_EnforceBounds.F90,v $
!   Revision 1.3  2008/12/06 08:44:40  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:52  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2004/01/29 22:59:29  haselbac
!   Initial revision
!
! ******************************************************************************







