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
! Purpose: Collection of routines related to Newton-Krylov schemes.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************

MODULE RFLU_ModPETScNewtonKrylov

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt,t_mixt_input
  USE ModMPI

  IMPLICIT NONE

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

  PRIVATE
  PUBLIC :: RFLU_PETSC_CreateVectors, & 
            RFLU_PETSC_CreateJacobian, &
            RFLU_PETSC_DestroyVectors, & 
            RFLU_PETSC_DestroyJacobian, &
            RFLU_PETSC_FormResidual, &
            RFLU_PETSC_FormResidualFirstOrder, &
            RFLU_PETSC_FormJacobian, &
            RFLU_PETSC_CreateColoring, &
            RFLU_PETSC_CreateApplicationOrdering, &
            RFLU_PETSC_DestroyApplicationOrdering, &
            RFLU_PETSC_GetDtScale

  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPETScNewtonKrylov.F90,v $ $Revision: 1.9 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  


! ******************************************************************************
!
! Purpose: Create and allocate memory for vectors needed by PETSc
!
! Description: None.
!
! Input: 
!   pRegion             Data for a grid region
!
! Output: 
!   x                   Solution vector (initial solution)
!   r                   Callback vector (zeros)
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_CreateVectors(pRegion)

  USE ModDataStruct, ONLY: t_region
  
  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"

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

  INTEGER :: ierr, ndof
  PetscScalar :: zero = 0.0
  TYPE(t_global),POINTER :: global  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_CreateVectors',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Count the number of degrees of freedom
! ******************************************************************************

  ndof = (CV_MIXT_ENER - CV_MIXT_DENS + 1) * pRegion%grid%nCells

! ******************************************************************************
! Create vectors
! ******************************************************************************

  CALL VecCreateMPIWithArray(PETSC_COMM_WORLD,ndof,PETSC_DECIDE, &
        pRegion%mixt%cv,pRegion%x,ierr)
  CALL VecDuplicate(pRegion%x,pRegion%r,ierr)
  CALL VecSet(pRegion%r,zero,ierr)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_CreateVectors




! ******************************************************************************
!
! Purpose: Create and allocate memory for matrices and solver needed by PETSc
!
! Description: None.
!
! Input: 
!   pRegion        data for grid region
!
! Output: 
!   A              the shell jacobian matrix for matrix-free operations
!   pA             the explicitely formed preconditioning jacobian matrix
!   snes           the PETSc SNES solver context
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_CreateJacobian(pRegion)

  USE ModDataStruct, ONLY: t_region
  
  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscsnes.h"

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

  INTEGER :: ierr, ndof, gndof, maxitsnes, maxitksp
  KSP :: ksp, subksp
  PC :: pc, subpc
  PetscScalar :: snestol, ksptol, snestrtol
  TYPE(t_global), POINTER :: global

  maxitsnes = 1
  maxitksp = 10000
  snestol = 1.0e-32   ! tight to avoid SNES thinking its converged
  ksptol = 1.0e-3     ! loose to facilitate fast solving
  snestrtol = 1.0e-32 ! tight to avoid SNES TR thinking its converged

! ******************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_CreateJacobian',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Create the Application Ordering context
! ******************************************************************************

  IF ( global%nProcAlloc > 1 ) THEN
    CALL RFLU_PETSC_CreateApplicationOrdering(pRegion)
  END IF ! global%nProcAlloc

! ******************************************************************************
! Count the number of degrees of freedom
! ******************************************************************************

  ndof = (CV_MIXT_ENER - CV_MIXT_DENS + 1) * pRegion%grid%nCells

! ******************************************************************************
! Set PETSc Trust Region Options
! ******************************************************************************

  CALL PetscOptionsSetValue('-snes_tr_delta0','1.e64',ierr)

! ******************************************************************************
! Find the global number of dof
! ******************************************************************************

  CALL MPI_ALLREDUCE(ndof,gndof,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! ******************************************************************************
! Create matrices
! ******************************************************************************

! ==============================================================================
! Create SNES solver context
! ==============================================================================

  CALL SNESCreate(PETSC_COMM_WORLD,pRegion%snes,ierr)

! ==============================================================================
! Create jacobian matrix and preconditioning matrix
! ==============================================================================

  CALL MatCreateSNESMF(pRegion%snes,pRegion%x,pRegion%A,ierr)
  CALL MatCreate(PETSC_COMM_WORLD,pRegion%preA,ierr)
  CALL MatSetSizes(pRegion%preA,ndof,ndof,gndof,gndof,ierr)
  CALL MatSetFromOptions(pRegion%preA,ierr)

! ==============================================================================
! Set solver parameters
! ==============================================================================

! ------------------------------------------------------------------------------
! Nonlinear solver parameters
! ------------------------------------------------------------------------------

  CALL SNESSetType(pRegion%snes,SNESTR,ierr)
  CALL SNESSetTolerances(pRegion%snes,snestol,PETSC_DEFAULT_DOUBLE_PRECISION, &
          PETSC_DEFAULT_DOUBLE_PRECISION,maxitsnes,PETSC_DEFAULT_INTEGER,ierr)
  CALL SNESSetTrustRegionTolerance(pRegion%snes,snestrtol,ierr)
  CALL SNESGetKSP(pRegion%snes,ksp,ierr)

! ------------------------------------------------------------------------------
! Linear solver & preconditioner parameters
! ------------------------------------------------------------------------------

  CALL KSPSetType(ksp,KSPGMRES,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCASM,ierr)
  CALL KSPSetTolerances(ksp,ksptol,PETSC_DEFAULT_DOUBLE_PRECISION, &
        PETSC_DEFAULT_DOUBLE_PRECISION,maxitksp,ierr)
  CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
  CALL KSPGMRESSetRestart(ksp,30,ierr)
  CALL KSPGMRESSetPreAllocateVectors(ksp,ierr)

! ------------------------------------------------------------------------------
! ASM Sub solver & preconditioner parameters
! ------------------------------------------------------------------------------

  CALL PetscOptionsSetValue('-sub_pc_type','ilu',ierr)
  CALL PetscOptionsSetValue('-sub_pc_ilu_levels','1',ierr)
  CALL PetscOptionsSetValue('-sub_ksp_type','preonly',ierr)

! ------------------------------------------------------------------------------
! Finalize solver parameters
! ------------------------------------------------------------------------------

  CALL SNESSetFromOptions(pRegion%snes,ierr)

! ==============================================================================
! Prefill Jacobian and Create Jacobian Coloring
! ==============================================================================

  IF ( global%myProcid == 0 ) THEN
     WRITE(*,*) 'Please wait while the Jacobian is prefilled and colored.  This may take a while...'
  END IF ! global%myProcid
  CALL RFLU_PETSC_CreateColoring(pRegion)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_CreateJacobian

  
! ******************************************************************************
!
! Purpose: Destroy the vectors that were created for PETSc
!
! Description: None.
!
! Input: 
!   pRegion           data for grid region
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_DestroyVectors(pRegion)

  USE ModDataStruct, ONLY: t_region  

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  INTEGER :: ierr
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_DestroyVectors',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Destroy the vectors
! ******************************************************************************

  CALL VecDestroy(pRegion%x,ierr)
  CALL VecDestroy(pRegion%r,ierr)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_DestroyVectors


! ******************************************************************************
!
! Purpose: Destroy the matrices and solver contexts that were created for PETSc
!
! Description: None.
!
! Input:
!   pRegion             data for grid region 
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_DestroyJacobian(pRegion)
  
  USE ModDataStruct, ONLY: t_region 

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  INTEGER :: ierr
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_DestroyJacobian',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Destroy the matrices and solver context
! ******************************************************************************
  
  IF ( global%nProcAlloc == 1 ) THEN
    CALL MatFDColoringDestroy(pRegion%fdcolor,ierr)
  END IF ! global%nProcAlloc
  CALL MatDestroy(pRegion%A,ierr)
  CALL MatDestroy(pRegion%preA,ierr)
  CALL SNESDestroy(pRegion%snes,ierr)

! ******************************************************************************
! Destroy the Application Ordering context
! ******************************************************************************

  IF ( global%nProcAlloc > 1 ) THEN
    CALL RFLU_PETSC_DestroyApplicationOrdering(pRegion)
  END IF ! global%nProcAlloc

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_DestroyJacobian






! ******************************************************************************
!
! Purpose: Explicitely forms the preconditioning Jacobian matrix for PETSc.
!
! Description: None.
!
! Input: 
!   snes            the PETSc SNES solver context
!   v               the vector at which to compute the jacobian
!   flag            PETSc flag denoting matrix structure
!   dummy           PETSc dummy context
!
! Output:
!   J               matrix-free jacobian matrix (is not modified)
!   pJ              explicitely formed preconditioning jacobian matrix
!   ierr            PETSc error code
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_FormJacobian(snes,v,J,pJ,flag,pRegion,ierr)

  USE ModDataStruct, ONLY: t_region
  
  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  SNES :: snes
  Mat :: J, pJ
  Vec :: v
  INTEGER :: ierr, order
  MatStructure :: flag
  TYPE(t_region),POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  REAL(kind=8) :: val
  INTEGER :: icg, ieq, neq, row
! DEBUG
  INTEGER :: n, i, k, low, high, nz, gnz
  PetscOffset xx_i
  PetscScalar xx_v(1)
! END DEBUG

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_FormJacobian',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Compute the preconditioning Jacobian matrix
! ******************************************************************************

  flag = DIFFERENT_NONZERO_PATTERN
  CALL SNESDefaultComputeJacobianColor(snes,v,J,pJ,flag, &
      pRegion%fdcolor,ierr)

! DEBUG
!
!  IF ( global%currentIter == 2 ) THEN
!   print*,'writing jacobian to file...'
!   if(global%myProcid == 0) then
!     OPEN(unit=123,file='jacobian_computational_0.m')
!   else
!    OPEN(unit=123,file='jacobian_computational_1.m')
!   endif
!   CALL VecGetSize(v,n,ierr)
!   CALL VecGetOwnershipRange(pRegion%x,low,high,ierr)
!   DO i = low, high-1
!     DO k = 1, n
!       CALL MatGetValues(pJ,1,i,1,k-1,val,ierr)
!       if(val /= 0.0) then
!          write(123,*) 'A(',i+1,',',k,') = ',val,';'
!       endif
!     ENDDO
!   ENDDO
!   close(123)
!   print*,'done writing jacobian to file'
!!   CALL VecView(v,PETSC_VIEWER_STDOUT_WORLD,ierr)
!   stop
!  END IF
!
!  IF ( global%myProcid == 0 ) THEN
!    DO i = 1, 5
!      DO k = 1, 5
!        CALL MatGetValues(pJ,1,i-1,1,k-1,val,ierr)
!        print*,'JAC',i,k,val
!      END DO
!    END DO
!  END IF
!  STOP
!
! END DEBUG

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_FormJacobian




! ******************************************************************************
!
! Purpose: Returns the residual of the nonlinear equations to PETSc.
!
! Description: None.
!
! Input: 
!   x                  vector at which to compute the residual
!   dummy              PETSc dummy input argument
!
! Output: 
!   f                  residual vector
!   ierr               PETSc error code
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_FormResidual(snes,x,f,pRegion,ierr)

  USE RFLU_ModResidual
  USE RFLU_ModMPI

  USE ModDataStruct, ONLY: t_region

  USE ModInterfaces, ONLY: RFLU_SetDependentVars, &
                           RFLU_CheckPositivity

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

  SNES snes
  Vec x, f
  INTEGER ierr, counter, icg, ieq, errorFlag
  PetscOffset xx_i, ff_i
  PetscScalar xx_v(1),ff_v(1),temp
  REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: cv, dv
  TYPE(t_region),POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_FormResidual',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ******************************************************************************
! Get the vectors
! ******************************************************************************

  CALL VecGetArray(x,xx_v,xx_i,ierr)
  CALL VecGetArray(f,ff_v,ff_i,ierr)

! ******************************************************************************
! Backup the original values for cv and dv
! ******************************************************************************

  ALLOCATE(cv(CV_MIXT_DENS:CV_MIXT_ENER,1:pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cv')
  END IF ! global%errorFlag

  ALLOCATE(dv(DV_MIXT_PRES:DV_MIXT_SOUN,1:pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dv')
  END IF ! global%errorFlag

  cv(:,:) = pRegion%mixt%cv(:,:)
  dv(:,:) = pRegion%mixt%dv(:,:)

! ******************************************************************************
! Compute copy the given vector into cv and check for validity
! ******************************************************************************

  counter = 0
  DO icg = 1, pRegion%grid%nCells
    DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
      counter = counter + 1
      pRegion%mixt%cv(ieq,icg) = xx_v(xx_i+counter)
    END DO ! ieq
  END DO ! icg

! ******************************************************************************
! Communicate virtual cell values & compute dependent variables
! ******************************************************************************

  CALL RFLU_MPI_ISendWrapper(pRegion)
!  CALL RFLU_MPI_CopyWrapper(regions)
  CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCells)
  CALL RFLU_MPI_RecvWrapper(pRegion)
  CALL RFLU_MPI_ClearRequestWrapper(pRegion)
  CALL RFLU_SetDependentVars(pRegion,pRegion%grid%nCells+1,pRegion%grid%nCellsTot)
  CALL RFLU_CheckPositivity(pRegion)

! ******************************************************************************
! Compute the residual
! ******************************************************************************

  CALL RFLU_RES_ComputeResidual(pRegion)

! ******************************************************************************
! Copy the residual array into the PETSc vector
! ******************************************************************************

  counter = 0

  DO icg = 1,pRegion%grid%nCells
    DO ieq = CV_MIXT_DENS,CV_MIXT_ENER
      counter = counter + 1
      IF ( global%flowType == FLOW_UNSTEADY ) THEN 
        temp = ( 1.5_RFREAL*pRegion%mixt%cv(ieq,icg) * pRegion%grid%vol(icg) &
               - 2.0_RFREAL*pRegion%mixt%cvOld1(ieq,icg) * pRegion%gridOld%vol(icg) &
               + 0.5_RFREAL*pRegion%mixt%cvOld2(ieq,icg) * pRegion%gridOld2%vol(icg) &
               ) / global%dtImposed
        ff_v(ff_i+counter) = pRegion%mixt%rhs(ieq,icg) &
          + pRegion%grid%vol(icg)/pRegion%dt(icg)*(pRegion%mixt%cv(ieq,icg)-cv(ieq,icg)) &
          + temp
      ELSE
        ff_v(ff_i+counter) = pRegion%mixt%rhs(ieq,icg) &
          + pRegion%grid%vol(icg)/pRegion%dt(icg)*(pRegion%mixt%cv(ieq,icg)-cv(ieq,icg))
      END IF ! global%flowType
    END DO ! ieq
  END DO ! icg

! ******************************************************************************
! Restore original values for cv and dv
! ******************************************************************************

  pRegion%mixt%cv(:,:) = cv(:,:)
  pRegion%mixt%dv(:,:) = dv(:,:)

  DEALLOCATE(cv,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cv')
  END IF !global%error

  DEALLOCATE(dv,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dv')
  END IF !global%error

! ******************************************************************************
! Put the vectors back together
! ******************************************************************************

  CALL VecRestoreArray(x,xx_v,xx_i,ierr)
  CALL VecRestoreArray(f,ff_v,ff_i,ierr)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_FormResidual








! ******************************************************************************
!
! Purpose: Returns the residual of the nonlinear equations to PETSc with first
!          order spacial accuracy.
!
! Description: None.
!
! Input: 
!   x                  vector at which to compute the residual
!   dummy              PETSc dummy input argument
!
! Output: 
!   f                  residual vector
!   ierr               PETSc error code
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_FormResidualFirstOrder(snes,x,f,pRegion,ierr)

  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

  SNES snes
  Vec x, f
  TYPE(t_region),POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  INTEGER :: ierr, oldSpaceOrder
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_FormResidualFirstOrder',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ==============================================================================
! Set the spacial order of accuracy to first order
! ==============================================================================

  oldSpaceOrder = pRegion%mixtInput%spaceOrder
  pRegion%mixtInput%spaceOrder = 1

! ==============================================================================
! Get the residual
! ==============================================================================

  CALL RFLU_PETSC_FormResidual(snes,x,f,pRegion,ierr)

! ==============================================================================
! Set the spacial order of accuracy to what it was to start with
! ==============================================================================

  pRegion%mixtInput%spaceOrder = oldSpaceOrder

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_FormResidualFirstOrder







! ******************************************************************************
!
! Purpose: Creates the Jacobian matrix coloring contexts needed by PETSc, and
!          preallocates and prefills the Jacobian matrix.
!
! Description: None.
!
! Input:
!   pRegion             data for grid region 
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_CreateColoring(pRegion)

  USE RFLU_ModColoring

  USE ModDataStruct, ONLY: t_region 
  USE RFLU_ModResidual, ONLY: RFLU_GetResidualSupport1, &
                              RFLU_GetResidualSupport2

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscsnes.h"
#include "finclude/petscis.h"
#include "finclude/petscao.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  INTEGER :: i, j, k, ifg, icg, ieq, irs, c1, c2, ierr, gndof, ndof, dof, color, neq
  INTEGER :: low, high, coloroffset, errorFlag, rsSizeMax, rsSize
  INTEGER,ALLOCATABLE,DIMENSION(:) :: nnzon, nnzoff, colors, pcolors, rs
  PetscScalar :: val
  ISColoring :: iscolor
  INTEGER,ALLOCATABLE,DIMENSION(:) :: rows, cols
  PetscScalar,ALLOCATABLE,DIMENSION(:,:) :: vals

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_CreateColoring',&
  'RFLU_ModPETScNewtonKrylov.F90')

! ==============================================================================
! Create the Connectivity Array
! ==============================================================================

  rsSizeMax = 1000

  ALLOCATE(rs(1:rsSizeMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'rs')
  END IF ! global%errorFlag

! ==============================================================================
! Preallocate the Jacobian Matrix
! ==============================================================================
 
! ------------------------------------------------------------------------------
! Get misc dimensions
! ------------------------------------------------------------------------------

  neq = CV_MIXT_ENER - CV_MIXT_DENS + 1
  CALL VecGetSize(pRegion%x,gndof,ierr)
  CALL VecGetLocalSize(pRegion%x,ndof,ierr)

! ------------------------------------------------------------------------------
! Allocate preallocation variables
! ------------------------------------------------------------------------------

  ALLOCATE(nnzon(1:ndof),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nnzon')
  END IF ! global%errorFlag

  ALLOCATE(nnzoff(1:ndof),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nnzoff')
  END IF ! global%errorFlag

! ------------------------------------------------------------------------------
! Calculate number of nonzeros in each row in the Jacobian. NOTE must be careful
! with Fortran and C-style numbering - nnzon uses Fortran numbering...
! ------------------------------------------------------------------------------

  nnzon(:) = neq
  nnzoff(:) = 0

  DO c1 = 1, pRegion%grid%nCells
! TEMPORARY
!    IF ( pRegion%mixtInput%spaceOrder == DISCR_ORDER_1 ) THEN
      CALL RFLU_GetResidualSupport1(pRegion,c1,rs,rsSizeMax,rsSize)
!    ELSE
!      CALL RFLU_GetResidualSupport2(pRegion,c1,rs,rsSizeMax,rsSize)
!    END IF ! pRegion%mixtInput%spaceOrder
! END TEMPORARY

    DO irs = 1, rsSize
      c2 = rs(irs)

      DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
        j = neq*(c1-1) + ieq 
        nnzon(j) = nnzon(j) + neq

        IF ( c2 <= pRegion%grid%nCells ) THEN
          j = neq*(c2-1) + ieq
          nnzon(j) = nnzon(j) + neq
        ELSE
          j = neq*(c1-1) + ieq 
          nnzoff(j) = nnzoff(j) + neq
        END IF ! c2
      END DO ! ieq
    END DO ! irs
  END DO ! c1

! ------------------------------------------------------------------------------
! Preallocate the matrix
! ------------------------------------------------------------------------------

  IF ( global%nProcAlloc == 1 ) THEN
    CALL MatSeqAIJSetPreallocation(pRegion%preA,0,nnzon,ierr)
  ELSE
    CALL MatMPIAIJSetPreallocation(pRegion%preA,0,nnzon,0,nnzoff,ierr)
  END IF ! global%nProcAlloc

! ------------------------------------------------------------------------------
! Deallocate preallocation variables
! ------------------------------------------------------------------------------

  DEALLOCATE(nnzon,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nnzon')
  END IF !global%error

  DEALLOCATE(nnzoff,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nnzoff')
  END IF !global%error

! ==============================================================================
! Pre-Fill the Jacobian Matrix
! ==============================================================================

! ------------------------------------------------------------------------------
! Allocate pre-fill variables
! ------------------------------------------------------------------------------

  ALLOCATE(vals(1:neq,1:neq),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vals')
  END IF ! global%errorFlag

  ALLOCATE(rows(1:neq),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'rows')
  END IF ! global%errorFlag

  ALLOCATE(cols(1:neq),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cols')
  END IF ! global%errorFlag

! ------------------------------------------------------------------------------
! Create the array of nonzero values
! ------------------------------------------------------------------------------

  DO i = 1, neq
    DO j = 1, neq
      vals(i,j) = 1.0
    END DO ! j
  END DO ! i

! ------------------------------------------------------------------------------
! Set nonzeros in off-diagonal blocks
! ------------------------------------------------------------------------------

  DO c1 = 1, pRegion%grid%nCells
! TEMPORARY
!    IF ( pRegion%mixtInput%spaceOrder == DISCR_ORDER_1 ) THEN
      CALL RFLU_GetResidualSupport1(pRegion,c1,rs,rsSizeMax,rsSize)
!    ELSE
!      CALL RFLU_GetResidualSupport2(pRegion,c1,rs,rsSizeMax,rsSize)
!    END IF ! pRegion%mixtInput%spaceOrder
! END TEMPORARY
    DO irs = 1, rsSize
      c2 = rs(irs)
      DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
        IF ( global%nProcAlloc == 1 ) THEN
          rows(ieq) = neq*(c1-1) + ieq - 1
          cols(ieq) = neq*(c2-1) + ieq - 1
        ELSE
          rows(ieq) = neq*(pRegion%grid%pc2sc(c1)-1) + ieq - 1
          cols(ieq) = neq*(pRegion%grid%pc2sc(c2)-1) + ieq - 1
        END IF ! global%nProcAlloc
      END DO ! ieq
      IF ( global%nProcAlloc > 1 ) THEN
        CALL AOApplicationToPetsc(pRegion%ao,neq,rows,ierr)
        CALL AOApplicationToPetsc(pRegion%ao,neq,cols,ierr)
      END IF ! global%nProcAlloc
      CALL MatSetValues(pRegion%preA,neq,rows,neq,cols,vals,INSERT_VALUES,ierr)
      IF ( c2 <= pRegion%grid%nCells ) THEN
        CALL MatSetValues(pRegion%preA,neq,cols,neq,rows,vals,INSERT_VALUES,ierr)
      END IF ! c2
    END DO ! irs
  END DO ! ifg

! ------------------------------------------------------------------------------
! Set nonzeros in on-diagonal blocks
! ------------------------------------------------------------------------------

  DO icg = 1, pRegion%grid%nCells
    DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
      IF ( global%nProcAlloc == 1 ) THEN
        rows(ieq) = neq*(icg-1) + ieq - 1
      ELSE
        rows(ieq) = neq*(pRegion%grid%pc2sc(icg)-1) + ieq - 1
      END IF ! global%nProcAlloc
    END DO ! ieq
    IF ( global%nProcAlloc > 1 ) THEN
      CALL AOApplicationToPetsc(pRegion%ao,neq,rows,ierr)
    END IF ! global%nProcAlloc
    CALL MatSetValues(pRegion%preA,neq,rows,neq,rows,vals,INSERT_VALUES,ierr)
  END DO ! icg  

! ------------------------------------------------------------------------------
! Deallocate pre-fill variables
! ------------------------------------------------------------------------------
  
  DEALLOCATE(vals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vals')
  END IF !global%error

  DEALLOCATE(rows,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'rows')
  END IF !global%error

  DEALLOCATE(cols,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cols')
  END IF !global%error

! ------------------------------------------------------------------------------
! Assemble the matrix
! ------------------------------------------------------------------------------

  CALL MatAssemblyBegin(pRegion%preA,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(pRegion%preA,MAT_FINAL_ASSEMBLY,ierr)

! DEBUG
!  print*,'writing jacobian to file...'
!  if(global%myProcid == 0) then
!    OPEN(unit=123,file='jacobian_computational_0.m')
!  else
!    OPEN(unit=123,file='jacobian_computational_1.m')
!  endif
!  CALL VecGetSize(pRegion%x,gndof,ierr)
!  CALL VecGetOwnershipRange(pRegion%x,low,high,ierr)
!  DO i = low, high-1
!    DO k = 1, gndof
!      CALL MatGetValues(pRegion%preA,1,i,1,k-1,val,ierr)
!      if(val /= 0.0) then
!         write(123,*) 'A(',i+1,',',k,') = ',val,';'
!      endif
!    ENDDO
!  ENDDO
!  close(123)
!  print*,'done writing jacobian to file'
!  stop
! END DEBUG

! ==============================================================================
! Create Coloring for Jacobian Evaluation
! ==============================================================================

! ------------------------------------------------------------------------------
! Allocate coloring variables
! ------------------------------------------------------------------------------

  CALL RFLU_COL_CreateColoring(pRegion)

  ALLOCATE(colors(1:ndof),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'colors')
  END IF ! global%errorFlag

! ------------------------------------------------------------------------------
! Read the cell coloring from the .col files
! ------------------------------------------------------------------------------

  CALL RFLU_COL_ReadColoring(pRegion)
  
! ------------------------------------------------------------------------------
! Compute the DOF coloring
! ------------------------------------------------------------------------------

  DO icg = 1, pRegion%grid%nCells
    DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
      colors(neq * (icg-1) + ieq) = neq * (pRegion%grid%col(icg) - 1) + ieq - 1
    END DO ! ieq
  END DO ! icg

! ------------------------------------------------------------------------------
! Give PETSc the coloring
! ------------------------------------------------------------------------------

  CALL ISColoringCreate(PETSC_COMM_WORLD,ndof,colors,iscolor,ierr)
  CALL MatFDColoringCreate(pRegion%preA,iscolor,pRegion%fdcolor,ierr)
  CALL MatFDColoringSetFromOptions(pRegion%fdcolor,ierr)
  CALL ISColoringDestroy(iscolor,ierr)

! ------------------------------------------------------------------------------
! Deallocate the coloring variables
! ------------------------------------------------------------------------------

  CALL RFLU_COL_DestroyColoring(pRegion)

  DEALLOCATE(colors,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'colors')
  END IF !global%error

! ==============================================================================
! Destroy the Connectivity Array
! ==============================================================================

  DEALLOCATE(rs,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'rs')
  END IF !global%error

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_CreateColoring




! ******************************************************************************
!
! Purpose: Creates the PETSc Application Ordering context to map Rocflu global
!          dof indices to PETSc global dof indices.
!
! Description: None.
!
! Input:
!   pRegion             data for grid region 
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_CreateApplicationOrdering(pRegion)
  
  USE ModDataStruct, ONLY: t_region 

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscao.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  INTEGER,ALLOCATABLE,DIMENSION(:) :: numpetsc, numapp
  INTEGER :: icg, ieq, neq, ipetsc, iapp, low, counter, ierr, errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_CreateApplicationOrdering',&
  'RFLU_ModPETScNewtonKrylov.F90')
  
! ==============================================================================
! Allocate the numbering arrays
! ============================================================================== 

  neq = CV_MIXT_ENER - CV_MIXT_DENS + 1  

  ALLOCATE(numpetsc(1:neq*pRegion%grid%nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'numpetsc')
  END IF ! global%errorFlag

  ALLOCATE(numapp(1:neq*pRegion%grid%nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'numapp')
  END IF ! global%errorFlag

! ==============================================================================
! Fill the application and PETSc numbering arrays
! ==============================================================================

  CALL VecGetOwnershipRange(pRegion%x,low,PETSC_NULL_INTEGER,ierr)
  ipetsc = low - 1
  counter = 0
  DO icg = 1, pRegion%grid%nCells
    DO ieq = CV_MIXT_DENS, CV_MIXT_ENER
      ipetsc = ipetsc + 1
      counter = counter + 1
      iapp = neq*(pRegion%grid%pc2sc(icg) - 1) + ieq - 1
      numapp(counter) = iapp
      numpetsc(counter) = ipetsc
    END DO ! ieq
  END DO ! icg

! ==============================================================================
! Create the Application Ordering context
! ==============================================================================

  CALL AOCreateBasic(PETSC_COMM_WORLD, neq*pRegion%grid%nCells, &
       numapp, numpetsc, pRegion%ao, ierr)

! ==============================================================================
! Deallocate the numbering arrays
! ============================================================================== 

  DEALLOCATE(numpetsc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'numpetsc')
  END IF !global%error

  DEALLOCATE(numapp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'numapp')
  END IF !global%error

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_CreateApplicationOrdering





! ******************************************************************************
!
! Purpose: Destroys the PETSc Application Ordering context
!
! Description: None.
!
! Input:
!   pRegion             data for grid region 
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_DestroyApplicationOrdering(pRegion)
  
  USE ModDataStruct, ONLY: t_region 

  IMPLICIT NONE

! ******************************************************************************
! Include PETSc headers
! ******************************************************************************

#include "finclude/petsc.h"
#include "finclude/petscao.h"

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  INTEGER :: ierr

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  CALL RegisterFunction(global,'RFLU_PETSC_DestroyApplicationOrdering',&
  'RFLU_ModPETScNewtonKrylov.F90')
  
! ==============================================================================
! Destroy the Application Ordering context
! ==============================================================================

  CALL AODestroy(pRegion%ao,ierr)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_DestroyApplicationOrdering






! ******************************************************************************
!
! Purpose: Provides the dt scaling factor for implicit runs
!
! Description: None.
!
! Input:
!   pRegion             data for grid region 
!
! Output: 
!   none
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_PETSC_GetDtScale(global,dtscale)
  
  USE ModDataStruct, ONLY: t_region 

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_global), POINTER :: global
  REAL(RFREAL) :: dtscale

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_PETSC_GetDtScale',&
  'RFLU_ModPETScNewtonKrylov.F90')
  
! ==============================================================================
! Compute the dt scaling
! ==============================================================================

  dtscale = 1.1 ** global%currentIter

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PETSC_GetDtScale






! ******************************************************************************
! END
! ******************************************************************************
  
END MODULE RFLU_ModPETScNewtonKrylov

















