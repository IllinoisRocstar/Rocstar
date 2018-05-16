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
! Purpose: Suite of routines to handle interaction with PETSc.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPETScAdmin.F90,v 1.4 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPETScAdmin

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC ::  RFLU_PETSC_Finalize, & 
             RFLU_PETSC_Init
      
#include "include/finclude/petsc.h"    
         
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPETScAdmin.F90,v $ $Revision: 1.4 $'               
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! ******************************************************************************
!
! Purpose: Finalize PETSc.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_Finalize(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global
    
! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_PETSC_Finalize',&
  'RFLU_ModPETScAdmin.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalizing PETSc...' 
    END IF ! global%verbLevel

! ******************************************************************************
!   Initialize PETSc
! ******************************************************************************

    CALL PetscFinalize(errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error    

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalizing PETSc done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_Finalize






! ******************************************************************************
!
! Purpose: Initialize PETSc.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_Init(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global
    
! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_PETSC_Init',&
  'RFLU_ModPETScAdmin.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing PETSc...' 
    END IF ! global%verbLevel

! ******************************************************************************
!   Initialize PETSc
! ******************************************************************************

    CALL PetscInitialize(PETSC_NULL_CHARACTER,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error    

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing PETSc done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_Init
  
  
  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModPETScAdmin


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPETScAdmin.F90,v $
! Revision 1.4  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/19 15:40:45  haselbac
! Initial revision
!
! ******************************************************************************








