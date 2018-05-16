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
! Purpose: Suite of routines to solve Poisson equation using PETSc.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPETScPoisson.F90,v 1.5 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPETScPoisson

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE RFLU_ModPETScAdmin

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_PETSC_BuildPoisson, & 
            RFLU_PETSC_CreatePoisson, & 
            RFLU_PETSC_SolvePressurePoisson, & 
            RFLU_PETSC_SetSolverPoisson   
 
#include "include/finclude/petsc.h"    
#include "include/finclude/petscvec.h"  
#include "include/finclude/petscmat.h"   
#include "include/finclude/petscksp.h"  
#include "include/finclude/petscpc.h"  
       
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPETScPoisson.F90,v $ $Revision: 1.5 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  





! ******************************************************************************
!
! Purpose: Build Poisson equation.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_BuildPoisson(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    KSP :: PETSC_ksp
    Mat :: PETSC_A
    MatNullSpace :: PETSC_nsp    
    PC :: PETSC_pc
    Vec :: PETSC_b,PETSC_x
    INTEGER :: PETSC_col,PETSC_row

    INTEGER :: c1,c2,errorFlag,i,icg,ifg,iPatch,j,jBeg,jEnd
    REAL(RFREAL) :: term
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PETSC_BuildPoisson',&
  'RFLU_ModPETScPoisson.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building Poisson...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ******************************************************************************
!   Initialize matrix
! ******************************************************************************
    
    DO i = 1,SIZE(pGrid%poissonA,1)
      pGrid%poissonA(i) = 0.0_RFREAL
    END DO ! i
  
    DO icg = 1,pGrid%nCells
      pGrid%poissonLast(icg) = pGrid%poissonFirst(icg) + 1
    END DO ! icg  
      
! ******************************************************************************
!   Assemble matrix
! ******************************************************************************

! ==============================================================================
!   Interior faces
! ==============================================================================
      
    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)
      
      term = pGrid%fn(4,ifg)/(pGrid%cofgDist(1,ifg) + pGrid%cofgDist(2,ifg))

      IF ( c1 <= pGrid%nCells ) THEN
        pGrid%poissonA(pGrid%poissonFirst(c1)) = & 
          pGrid%poissonA(pGrid%poissonFirst(c1)) + term
        pGrid%poissonA(pGrid%poissonLast(c1)) = -term
        pGrid%poissonLast(c1) = pGrid%poissonLast(c1) + 1
      END IF ! c1

      IF ( c2 <= pGrid%nCells ) THEN
        pGrid%poissonA(pGrid%poissonFirst(c2)) = & 
          pGrid%poissonA(pGrid%poissonFirst(c2)) + term
        pGrid%poissonA(pGrid%poissonLast(c2)) = -term
        pGrid%poissonLast(c2) = pGrid%poissonLast(c2) + 1
      END IF ! c2
    END DO ! ifg

! ==============================================================================
!   Boundary faces
! ==============================================================================

! TO BE DONE    
!   loop over Dirichlet boundary patches
!   but generally we do not have Dirichlet boundary conditions for pressure
! END TO BE DONE

! ******************************************************************************
!   Set matrix. NOTE matrix lower boundary are always zero, even in FORTRAN.
! ******************************************************************************

    PETSC_A = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_A)

    DO icg = 1,pGrid%nCells
      jBeg = pGrid%poissonFirst(icg)
      jEnd = pGrid%poissonFirst(icg) + pGrid%poissonNNZ(icg)-1
      
      DO j = jBeg,jEnd
        PETSC_row = icg - 1
        PETSC_col = pGrid%poissonCol(j) - 1


        CALL MatSetValues(PETSC_A,1,PETSC_row,1,PETSC_col, &
                          pGrid%poissonA(j),INSERT_VALUES,errorFlag)

        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
        END IF ! global%error
      END DO ! j 
    END DO ! icg

    CALL MatAssemblyBegin(PETSC_A,MAT_FINAL_ASSEMBLY,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
    
    CALL MatAssemblyEnd(PETSC_A,MAT_FINAL_ASSEMBLY,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error        
  
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building Poisson done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_BuildPoisson







! ******************************************************************************
!
! Purpose: 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_CreatePoisson(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    Mat :: PETSC_A
    Vec :: PETSC_b,PETSC_x

    INTEGER :: c1,c2,errorFlag,icg,ifg,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PETSC_CreatePoisson',&
  'RFLU_ModPETScPoisson.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Poisson...' 
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ******************************************************************************
!   Allocate memory for non-zero pattern and access arrays
! ******************************************************************************

    ALLOCATE(pGrid%poissonNNZ(pGrid%nCells),STAT=errorFlag) 
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%poissonNNZ')
    END IF ! global%error 

    ALLOCATE(pGrid%poissonFirst(pGrid%nCells),STAT=errorFlag) 
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%poissonFirst')
    END IF ! global%error 

    ALLOCATE(pGrid%poissonLast(pGrid%nCells),STAT=errorFlag) 
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%poissonLast')
    END IF ! global%error     
    
! ******************************************************************************
!   Set non-zero pattern and access arrays
! ******************************************************************************
    
    DO icg = 1,pGrid%nCells
      pGrid%poissonNNZ(icg) = 1
    END DO ! icg
    
    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      IF ( c1 <= pGrid%nCells ) THEN 
        pGrid%poissonNNZ(c1) = pGrid%poissonNNZ(c1) + 1
      END IF ! c1
      
      IF ( c2 <= pGrid%nCells ) THEN 
        pGrid%poissonNNZ(c2) = pGrid%poissonNNZ(c2) + 1
      END IF ! c2
    END DO ! ifg

    pGrid%poissonFirst(1) = 1

    DO icg = 2,pGrid%nCells
      pGrid%poissonFirst(icg) = pGrid%poissonFirst(icg-1) & 
                              + pGrid%poissonNNZ(icg-1)
    END DO ! icg

    ALLOCATE(pGrid%poissonA(SUM(pGrid%poissonNNZ)),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%poissonA')
    END IF ! global%error
    
    ALLOCATE(pGrid%poissonCol(SUM(pGrid%poissonNNZ)),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%poissonCol')
    END IF ! global%error    
    
    DO icg = 1,pGrid%nCells
      pGrid%poissonCol(pGrid%poissonFirst(icg)) = icg
      
      pGrid%poissonLast(icg) = pGrid%poissonFirst(icg) + 1
    END DO ! icg 
    
    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      IF ( c1 <= pGrid%nCells ) THEN
        pGrid%poissonCol(pGrid%poissonLast(c1)) = c2
        pGrid%poissonLast(c1) = pGrid%poissonLast(c1) + 1
      END IF ! c1
      
      IF (c2 <= pGrid%nCells) THEN
        pGrid%poissonCol(pGrid%poissonLast(c2)) = c1
        pGrid%poissonLast(c2) = pGrid%poissonLast(c2) + 1
      END IF ! c2
    END DO ! ifg

! ******************************************************************************
!   Set matrix structure in PETSc
! ******************************************************************************

    CALL MatCreateSeqAIJ(PETSC_COMM_SELF,pGrid%nCells,pGrid%nCells, &
                         PETSC_NULL_INTEGER,pGrid%poissonNNZ,PETSC_A, & 
                         errorFlag)

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error                             

    CALL VecCreateSeqWithArray(PETSC_COMM_SELF,pGrid%nCells,PETSC_NULL_SCALAR, &
                               PETSC_b,errorFlag)

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 

    CALL VecDuplicate(PETSC_b,PETSC_x,errorFlag)

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 

    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_A) = PETSC_A
    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_B) = PETSC_b
    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_X) = PETSC_x            

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Poisson done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_CreatePoisson
  
  
  
  



  
! ******************************************************************************
!
! Purpose: Set solver options for Poisson equation.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_SetSolverPoisson(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    KSP :: PETSC_ksp
    Mat :: PETSC_A
    PC :: PETSC_pc
    MatNullSpace :: PETSC_nsp

    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PETSC_SetSolverPoisson',&
  'RFLU_ModPETScPoisson.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting solver Poisson...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set solver and preconditioner options
! ******************************************************************************
       
    PETSC_A = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_A)       
       
    CALL KSPCreate(PETSC_COMM_SELF,PETSC_ksp,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
    
    CALL KSPSetOperators(PETSC_ksp,PETSC_A,PETSC_A,SAME_NONZERO_PATTERN, & 
                         errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 

    CALL KSPSetFromOptions(PETSC_ksp,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
       
! TEMPORARY, needs to be sorted out 
!    CALL MatNullSpaceCreate(PETSC_COMM_SELF,PETSC_TRUE,0,PETSC_NULL_VEC, & 
!                            PETSC_nsp,errorFlag)
! END TEMPORARY

! Addition by hzhao
    CALL MatNullSpaceCreate(PETSC_COMM_SELF,PETSC_TRUE,0,0, & 
                            PETSC_nsp,errorFlag)
! End addition by hzhao

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 

    CALL KSPSetNullSpace(PETSC_ksp,PETSC_nsp,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 

    CALL KSPSetType(PETSC_ksp,KSPGMRES,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
        
    CALL KSPGetPC(PETSC_ksp,PETSC_pc,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
        
    CALL PCSetType(PETSC_pc,PCILU,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error                                                
  
    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_KSP) = PETSC_ksp
    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_PC ) = PETSC_pc
    pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_NSP) = PETSC_nsp                 
  
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting solver Poisson done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_SetSolverPoisson
  
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Solve Poisson equation.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PETSC_SolvePressurePoisson(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    KSP :: PETSC_ksp
    Mat :: PETSC_A
    PC :: PETSC_pc
    Vec :: PETSC_b,PETSC_x    

    INTEGER :: c1,c2,errorFlag,ifg,ifl,iPatch
    REAL(RFREAL) :: flx
    PetscReal, DIMENSION(:), ALLOCATABLE :: b,x
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
! Addition by hzhao
    INTEGER :: cvMixtPres
! End addition by hzhao
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PETSC_SolvePressurePoisson',&
  'RFLU_ModPETScPoisson.F90')

    pGrid => pRegion%grid  

! ******************************************************************************
!   Allocate memory
! ******************************************************************************
 
    ALLOCATE(x(pGrid%nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x')
    END IF ! global%error

    ALLOCATE(b(pGrid%nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'b')
    END IF ! global%error

! ******************************************************************************
!   Set right-hand side
! ******************************************************************************              
    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

! Comment by hzhao
!      flx = pRegion%mixt%vfMixt(ifg)*pGrid%fn(XCOORD,ifg) & 
!          + pRegion%mixt%vfMixt(ifg)*pGrid%fn(YCOORD,ifg) & 
!          + pRegion%mixt%vfMixt(ifg)*pGrid%fn(ZCOORD,ifg)
! End comment by hzhao
! Addition by hzhao
      flx = pRegion%mixt%vfMixt(ifg)*pGrid%fn(4,ifg)
! End addition by hzhao
          
      IF ( c1 <= pGrid%nCells ) THEN 
        b(c1) = b(c1) + flx
      END IF ! c1
      
      IF ( c2 <= pGrid%nCells ) THEN 
        b(c2) = b(c2) - flx
      END IF ! c2                 
    END DO ! ifg

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)            
      
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

! Comment by hzhao        
!        flx = pPatch%vfMixt(ifl)*pPatch%fn(XCOORD,ifl) & 
!            + pPatch%vfMixt(ifl)*pPatch%fn(YCOORD,ifl) & 
!            + pPatch%vfMixt(ifl)*pPatch%fn(ZCOORD,ifl)
! End comment by hzhao        
! Addition by hzhao
        flx = pPatch%vfMixt(ifl)*pPatch%fn(4,ifl) & 
! End addition by hzhao
            
        b(c1) = b(c1) + flx
      END DO ! ifl
    END DO ! iPatch

! ******************************************************************************
!   Set left-hand side
! ******************************************************************************              

! TO DO 
! 
! END TO DO 

! ******************************************************************************
!   Solve
! ******************************************************************************              

    PETSC_A = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_A)
    PETSC_b = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_B)
    PETSC_x = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_X)  
    
    PETSC_ksp = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_KSP)
    PETSC_pc  = pRegion%poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_PC )
     
    CALL VecPlaceArray(PETSC_b,b,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error 
    
    CALL VecPlaceArray(PETSC_x,x,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error      

! Addition by hzhao
! Solvability condition
    b = b - SUM(b)/SIZE(b,1)
! End addition by hzhao

    CALL KSPSolve(PETSC_ksp,PETSC_b,PETSC_x,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error    

! ******************************************************************************
!   Store pressure
! ******************************************************************************              
! Addition by hzhao
    cvMixtPres = RFLU_GetCvLoc(global, FLUID_MODEL_INCOMP, CV_MIXT_PRES)
    pRegion%mixt%cv(cvMixtPres,1:pGrid%nCells) = x
! End addition by hzhao

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
 
    DEALLOCATE(x,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x')
    END IF ! global%error

    DEALLOCATE(b,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'b')
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_SolvePressurePoisson
  
  
  
  
  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModPETScPoisson


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPETScPoisson.F90,v $
! Revision 1.5  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.2  2005/01/13 21:39:38  haselbac
! Incorporated Hongs additions for testing
!
! Revision 1.1  2004/12/19 15:40:45  haselbac
! Initial revision
!
! ******************************************************************************










