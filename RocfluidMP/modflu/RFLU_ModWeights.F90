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
! Purpose: Suite of routines to compute stencil weights.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModWeights.F90,v 1.16 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModWeights

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI
    
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_CreateWtsBF2CWrapper, &
            RFLU_CreateWtsC2CWrapper, &
            RFLU_CreateWtsF2CWrapper, &
            RFLU_ComputeWtsBF2CWrapper, &
            RFLU_ComputeWtsC2CWrapper, &
            RFLU_ComputeWtsF2CWrapper, &
            RFLU_ComputeWtsX2C_1D, &
            RFLU_DestroyWtsBF2CWrapper, & 
            RFLU_DestroyWtsC2CWrapper, &
            RFLU_DestroyWtsF2CWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModWeights.F90,v $ $Revision: 1.16 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  



! *******************************************************************************
!
! Purpose: Create weights for boundary face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!   orderInput  Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_CreateWtsBF2C(pRegion,pPatch,orderInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderInput  
    TYPE(t_patch), POINTER :: pPatch         
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl,order
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid      

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsBF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating boundary face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Modify order so that can run with first-order scheme
! ******************************************************************************

    order = MAX(orderInput,1)

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsBF2C(pRegion,pPatch)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    SELECT CASE ( order ) 
    
! ==============================================================================
!     Linear approximation
! ==============================================================================    
    
      CASE ( 1 ) 
        IF ( pPatch%bcType /= BC_VIRTUAL ) THEN  
          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 1 ) ! NOTE should never reach here
            CASE ( 2 ) 
              DO ifl = 1,pPatch%nBFaces
                ALLOCATE(pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11:XYZ_MOM_33), & 
                         STAT=errorFlag)  
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                                 'pPatch%bf2cs%xyzMoms')
                END IF ! global%error
              END DO ! ifl            
            CASE ( 3 )            
              DO ifl = 1,pPatch%nBFaces
                ALLOCATE(pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11:XYZ_MOM_44), & 
                         STAT=errorFlag)  
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                                 'pPatch%bf2cs%xyzMoms')
                END IF ! global%error
              END DO ! ifl
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! dimens
        END IF ! pPatch%bcType            

! ==============================================================================
!     Default
! ==============================================================================    

      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating boundary face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsBF2C  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Wrapper routine for creating boundary-face-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch 	Pointer to patch
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateWtsBF2CWrapper(pRegion,pPatch,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsBF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensBFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_CreateWtsBF2C(pRegion,pPatch,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsBF2CWrapper
  
  
  
  
  


! *******************************************************************************
!
! Purpose: Create weights for cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_CreateWtsC2C(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order  
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsC2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating cell-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsC2C(pRegion)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    SELECT CASE ( order ) 

! ==============================================================================
!     Linear approximation
! ==============================================================================    

      CASE ( 1 )
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) 
          CASE ( 2 )
            DO icg = 1,pGrid%nCellsTot 
              ALLOCATE(pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11:XYZ_MOM_33), & 
                       STAT=errorFlag)                   
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                               'pGrid%c2cs%xyzMoms')
              END IF ! global%error                  
            END DO ! icg           
          CASE ( 3 ) 
            DO icg = 1,pGrid%nCellsTot 
              ALLOCATE(pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11:XYZ_MOM_44), & 
                       STAT=errorFlag)                   
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                               'pGrid%c2cs%xyzMoms')
              END IF ! global%error                  
            END DO ! icg
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%dimens  

! ==============================================================================
!     Default
! ==============================================================================    

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating cell-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsC2C








! *******************************************************************************
!
! Purpose: Wrapper routine for creating cell-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateWtsC2CWrapper(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsC2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensCells )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_CreateWtsC2C(pRegion,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsC2CWrapper








! *******************************************************************************
!
! Purpose: Create weights for face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   orderInput  Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_CreateWtsF2C(pRegion,orderInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderInput  
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifg,order
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid            

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Modify order so that can run with first-order scheme
! ******************************************************************************

    order = MAX(orderInput,1)

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsF2C(pRegion)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    SELECT CASE ( order ) 
    
! ==============================================================================
!     Linear approximation
! ==============================================================================    
    
      CASE ( 1 ) 
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) ! NOTE should never reach here
          CASE ( 2 ) 
            DO ifg = 1,pGrid%nFaces
              ALLOCATE(pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11:XYZ_MOM_33), & 
                       STAT=errorFlag)                   
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                               'pGrid%f2cs%xyzMoms')
              END IF ! global%error        
            END DO ! ifg          
          CASE ( 3 )           
            DO ifg = 1,pGrid%nFaces
              ALLOCATE(pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11:XYZ_MOM_44), & 
                       STAT=errorFlag)                   
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                               'pGrid%f2cs%xyzMoms')
              END IF ! global%error        
            END DO ! ifg
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! dimens
    
! ==============================================================================
!     Default
! ==============================================================================    

      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsF2C







! *******************************************************************************
!
! Purpose: Wrapper routine for creating face-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateWtsF2CWrapper(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateWtsF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_CreateWtsF2C(pRegion,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateWtsF2CWrapper






! *******************************************************************************
!
! Purpose: Compute weights for boundary face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!   orderInput  Desired order 
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ComputeWtsBF2C(pRegion,pPatch,orderInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderInput  
    TYPE(t_patch), POINTER :: pPatch    
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,nMembs,icg,ifl,isl,order
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: dr    
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid        

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsBF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Computing boundary face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Modify order so that can run with first-order scheme
! ******************************************************************************

    order = MAX(orderInput,1)
    
! ******************************************************************************
!   Compute weights
! ******************************************************************************

    SELECT CASE ( order ) 
    
! ==============================================================================
!     Linear approximation
! ==============================================================================    
    
      CASE ( 1 )         
        IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 1 ) 
            CASE ( 2 ) 
              DO ifl = 1,pPatch%nBFaces
                nMembs = pPatch%bf2cs(ifl)%nCellMembs

                ALLOCATE(dr(XCOORD:YCOORD,nMembs),STAT=errorFlag)
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
                END IF ! global%error

                DO isl = 1,nMembs
                  icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                  dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg) & 
                                 - pPatch%fc(XCOORD,ifl)
                  dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg) & 
                                 - pPatch%fc(YCOORD,ifl)
                END DO ! isl        

                CALL RFLU_ComputeStencilMoments2D1(global,nMembs,dr, & 
                                                   pPatch%bf2cs(ifl)%xyzMoms)

                DEALLOCATE(dr,STAT=errorFlag)
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
                END IF ! global%error              
              END DO ! ifl            
            CASE ( 3 ) 
              DO ifl = 1,pPatch%nBFaces
                nMembs = pPatch%bf2cs(ifl)%nCellMembs

                ALLOCATE(dr(XCOORD:ZCOORD,nMembs),STAT=errorFlag)
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
                END IF ! global%error

                DO isl = 1,nMembs
                  icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                  dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg) & 
                                 - pPatch%fc(XCOORD,ifl)
                  dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg) & 
                                 - pPatch%fc(YCOORD,ifl)
                  dr(ZCOORD,isl) = pGrid%cofg(ZCOORD,icg) & 
                                 - pPatch%fc(ZCOORD,ifl)                    
                END DO ! isl        

                CALL RFLU_ComputeStencilMoments3D1(global,nMembs,dr, & 
                                                   pPatch%bf2cs(ifl)%xyzMoms)

                DEALLOCATE(dr,STAT=errorFlag)
                global%error = errorFlag
                IF ( global%error /= ERR_NONE ) THEN 
                  CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
                END IF ! global%error              
              END DO ! ifl
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%dimens
        END IF ! pPatch%bcType 

! ==============================================================================
!     Default
! ==============================================================================    
        
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Computing boundary face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsBF2C








! *******************************************************************************
!
! Purpose: Wrapper routine for computing boundary-face-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch 
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsBF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensBFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_ComputeWtsBF2C(pRegion,pPatch,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsBF2CWrapper







! *******************************************************************************
!
! Purpose: Compute weights for cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ComputeWtsC2C(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order  
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,nMembs,icg,icg2,isl
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: dr    
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsC2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Computing cell-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute weights
! ******************************************************************************

    SELECT CASE ( order ) 

! ==============================================================================
!     Linear approximation
! ==============================================================================    

      CASE ( 1 ) 
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) ! NOTE should never reach here
          CASE ( 2 ) 
            DO icg = 1,pGrid%nCellsTot
              nMembs = pGrid%c2cs(icg)%nCellMembs
            
              ALLOCATE(dr(XCOORD:YCOORD,nMembs),STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
              END IF ! global%error

              DO isl = 1,nMembs
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg2)-pGrid%cofg(XCOORD,icg)
                dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg2)-pGrid%cofg(YCOORD,icg)
              END DO ! isl        

              CALL RFLU_ComputeStencilMoments2D1(global,nMembs,dr, & 
                                                 pGrid%c2cs(icg)%xyzMoms)

              DEALLOCATE(dr,STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
              END IF ! global%error        
            END DO ! icg          
          CASE ( 3 )           
            DO icg = 1,pGrid%nCellsTot
              nMembs = pGrid%c2cs(icg)%nCellMembs            
            
              ALLOCATE(dr(XCOORD:ZCOORD,nMembs),STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
              END IF ! global%error

              DO isl = 1,nMembs
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg2)-pGrid%cofg(XCOORD,icg)
                dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg2)-pGrid%cofg(YCOORD,icg)
                dr(ZCOORD,isl) = pGrid%cofg(ZCOORD,icg2)-pGrid%cofg(ZCOORD,icg)                    
              END DO ! isl        

              CALL RFLU_ComputeStencilMoments3D1(global,nMembs,dr, & 
                                                 pGrid%c2cs(icg)%xyzMoms)

              DEALLOCATE(dr,STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
              END IF ! global%error        
            END DO ! icg
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%dimens

! ==============================================================================
!     Default
! ==============================================================================    

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Computing cell-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsC2C








! *******************************************************************************
!
! Purpose: Wrapper routine for computing cell-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeWtsC2CWrapper(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsC2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensCells )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_ComputeWtsC2C(pRegion,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsC2CWrapper







! *******************************************************************************
!
! Purpose: Compute weights for face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   orderInput  Desired order 
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ComputeWtsF2C(pRegion,orderInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderInput  
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,nMembs,icg,ifg,isl,order
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: dr    
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Computing face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Modify order so that can run with first-order scheme
! ******************************************************************************

    order = MAX(orderInput,1)
    
! ******************************************************************************
!   Compute weights
! ******************************************************************************

    SELECT CASE ( order ) 
    
! ==============================================================================
!     Linear approximation
! ==============================================================================    
    
      CASE ( 1 ) 
      
! ------------------------------------------------------------------------------
!       Interior faces
! ------------------------------------------------------------------------------      
      
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) ! NOTE should never reach here
          CASE ( 2 )
            DO ifg = 1,pGrid%nFaces
              nMembs = pGrid%f2cs(ifg)%nCellMembs
              
              ALLOCATE(dr(XCOORD:YCOORD,nMembs),STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
              END IF ! global%error

              DO isl = 1,nMembs
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
                dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)
              END DO ! isl        

              CALL RFLU_ComputeStencilMoments2D1(global,nMembs,dr, & 
                                                 pGrid%f2cs(ifg)%xyzMoms)

              DEALLOCATE(dr,STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
              END IF ! global%error        
            END DO ! ifg           
          CASE ( 3 ) 
            DO ifg = 1,pGrid%nFaces
              nMembs = pGrid%f2cs(ifg)%nCellMembs            
            
              ALLOCATE(dr(XCOORD:ZCOORD,nMembs),STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
              END IF ! global%error

              DO isl = 1,nMembs
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg)-pGrid%fc(XCOORD,ifg)
                dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg)-pGrid%fc(YCOORD,ifg)
                dr(ZCOORD,isl) = pGrid%cofg(ZCOORD,icg)-pGrid%fc(ZCOORD,ifg)                    
              END DO ! isl        

              CALL RFLU_ComputeStencilMoments3D1(global,nMembs, & 
                                                 dr,pGrid%f2cs(ifg)%xyzMoms)

              DEALLOCATE(dr,STAT=errorFlag)
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) THEN 
                CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
              END IF ! global%error        
            END DO ! ifg
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%dimens

! ==============================================================================
!     Default
! ==============================================================================    
        
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! order

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Computing face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsF2C







! *******************************************************************************
!
! Purpose: Wrapper routine for computing face-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   order       Desired order
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeWtsF2CWrapper(pRegion,order)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: order
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeWtsF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_ComputeWtsF2C(pRegion,order)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsF2CWrapper







! *******************************************************************************
!
! Purpose: Compute weights for 1D x-to-cell stencil.
!
! Description: Compute Lagrangian polynomial stencil weights from Fornberg 
!   algorithm.
!
! Input:
!   global      Pointer to global data
!   m           Order of derivative for which weights are sought
!   nMembs      Number of stencil members
!   x           Locations of stencil members
!   z           Location at which weights are sough
!
! Output:
!   w           Stencil weights (corresponding to c array in Fornberg)
!
! Notes: 
!   1. See Fornberg, SIAM Rev., Vol.40, No.3, pp.685-691, Sep. 1998.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_ComputeWtsX2C_1D(global,m,nMembs,x,z,w) 

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: m,nMembs
    REAL(RFREAL), INTENT(IN) :: z
    REAL(RFREAL), INTENT(IN) :: x(0:nMembs-1)
    REAL(RFREAL), INTENT(OUT) :: w(0:nMembs-1)
    TYPE(t_global), POINTER :: global  

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,i,j,k,mn,n
    REAL(RFREAL) :: c1,c2,c3,c4,c5 ! Same notation as Fornberg
    REAL(RFREAL) :: c(0:nMembs-1,0:m)    

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeWtsX2C_1D',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    n = nMembs-1

! TO DO 
!   Add check for order (given nMembs)
! END TO DO  

! ******************************************************************************
!   Initialize
! ******************************************************************************
    
    c1 = 1.0_RFREAL
    c4 = x(0) - z

    DO k = 0,m
      DO j = 0,n
        c(j,k) = 0.0_RFREAL
      END DO ! j
    END DO ! k

    c(0,0) = 1.0_RFREAL

! ******************************************************************************
!   Compute weights
! ******************************************************************************

    DO i = 1,n
      mn = MIN(i,m)
      c2 = 1.0_RFREAL
      c5 = c4
      c4 = x(i) - z
      
      DO j = 0,i-1
        c3 = x(i) - x(j)
        c2 = c2*c3
        
        IF ( j == (i-1) ) THEN
          DO k = mn,1,-1
            c(i,k) = c1*(k*c(i-1,k-1) - c5*c(i-1,k))/c2
          END DO ! k
          
          c(i,0) = -c1*c5*c(i-1,0)/c2
        END IF ! j
        
        DO k = mn,1,-1
          c(j,k) = (c4*c(j,k) - k*c(j,k-1))/c3
        END DO ! k
        
        c(j,0) = c4*c(j,0)/c3
      END DO ! j
      
      c1 = c2
    END DO ! i

    DO j = 0,nMembs-1
      w(j) = c(j,m)
    END DO ! j
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_ComputeWtsX2C_1D






! *******************************************************************************
!
! Purpose: Compute weights from coordinate differences.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   nMembs      Number of members in stencil
!   dr          Coordinate differences
!
! Output: 
!   xyzMoms     Coordinate moments
!
! Notes:
!   1. Use inverse-distance weighting.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeStencilMoments2D1(global,nMembs,dr,xyzMoms)
  
    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: nMembs
    REAL(RFREAL), INTENT(INOUT) :: dr(XCOORD:YCOORD,nMembs)
    REAL(RFREAL), INTENT(INOUT) :: xyzMoms(XYZ_MOM_11:XYZ_MOM_33)
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: isl  
    REAL(RFREAL) :: dx,dy,ir11,ir22,ir33,r11,r12,r13,r22,r23, &
                    r33,wt
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeStencilMoments2D1',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Compute moments
! ******************************************************************************
                
! ==============================================================================
!   Initialize weights
! ==============================================================================    
    
    xyzMoms(XYZ_MOM_11) = 0.0_RFREAL

    xyzMoms(XYZ_MOM_12) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_22) = 0.0_RFREAL

    xyzMoms(XYZ_MOM_13) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_23) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_33) = 0.0_RFREAL

! ==============================================================================
!   Compute weights
! ==============================================================================

    DO isl = 1,nMembs          
      dx = dr(XCOORD,isl) 
      dy = dr(YCOORD,isl)                      

      wt = 1.0_RFREAL/SQRT(dx*dx + dy*dy)

      dx = wt*dx 
      dy = wt*dy

      xyzMoms(XYZ_MOM_11) = xyzMoms(XYZ_MOM_11) + dx*dx

      xyzMoms(XYZ_MOM_12) = xyzMoms(XYZ_MOM_12) + dx*dy
      xyzMoms(XYZ_MOM_22) = xyzMoms(XYZ_MOM_22) + dy*dy

      xyzMoms(XYZ_MOM_13) = xyzMoms(XYZ_MOM_13) + wt*dx
      xyzMoms(XYZ_MOM_23) = xyzMoms(XYZ_MOM_23) + wt*dy
      xyzMoms(XYZ_MOM_33) = xyzMoms(XYZ_MOM_33) + wt*wt  
    END DO ! isl

    r11  = SQRT(xyzMoms(XYZ_MOM_11)) 
    ir11 = 1.0_RFREAL/r11

    r12  = ir11*xyzMoms(XYZ_MOM_12)
    r22  = SQRT(xyzMoms(XYZ_MOM_22) - r12*r12)
    ir22 = 1.0_RFREAL/r22

    r13  = ir11*xyzMoms(XYZ_MOM_13)
    r23  = ir22*(xyzMoms(XYZ_MOM_23) -  r12*r13           ) 
    r33  = SQRT(xyzMoms(XYZ_MOM_33)  - (r13*r13 + r23*r23))
    ir33 = 1.0_RFREAL/r33

! ==============================================================================
!   Store weights
! ==============================================================================

    xyzMoms(XYZ_MOM_11) = r11        

    xyzMoms(XYZ_MOM_12) = r12
    xyzMoms(XYZ_MOM_22) = r22

    xyzMoms(XYZ_MOM_13) = r13
    xyzMoms(XYZ_MOM_23) = r23
    xyzMoms(XYZ_MOM_33) = r33 
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)    
  
  END SUBROUTINE RFLU_ComputeStencilMoments2D1








! *******************************************************************************
!
! Purpose: Compute weights from coordinate differences.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   nMembs      Number of members in stencil
!   dr          Coordinate differences
!
! Output: 
!   xyzMoms     Coordinate moments
!
! Notes: 
!   1. Use inverse-distance weighting.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeStencilMoments3D1(global,nMembs,dr,xyzMoms)
  
    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: nMembs
    REAL(RFREAL), INTENT(INOUT) :: dr(XCOORD:ZCOORD,nMembs)
    REAL(RFREAL), INTENT(INOUT) :: xyzMoms(XYZ_MOM_11:XYZ_MOM_44)
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: isl  
    REAL(RFREAL) :: dx,dy,dz,ir11,ir22,ir33,r11,r12,r13,r14,r22,r23,r24, &
                    r33,r34,r44,wt
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeStencilMoments3D1',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Compute moments
! ******************************************************************************
                
! ==============================================================================
!   Initialize weights
! ==============================================================================    
    
    xyzMoms(XYZ_MOM_11) = 0.0_RFREAL

    xyzMoms(XYZ_MOM_12) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_22) = 0.0_RFREAL

    xyzMoms(XYZ_MOM_13) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_23) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_33) = 0.0_RFREAL

    xyzMoms(XYZ_MOM_14) = 0.0_RFREAL
    xyzMoms(XYZ_MOM_24) = 0.0_RFREAL                              
    xyzMoms(XYZ_MOM_34) = 0.0_RFREAL 
    xyzMoms(XYZ_MOM_44) = 0.0_RFREAL                                             

! ==============================================================================
!   Compute weights
! ==============================================================================

    DO isl = 1,nMembs          
      dx = dr(XCOORD,isl) 
      dy = dr(YCOORD,isl) 
      dz = dr(ZCOORD,isl)                     

      wt = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)

      dx = wt*dx 
      dy = wt*dy
      dz = wt*dz

      xyzMoms(XYZ_MOM_11) = xyzMoms(XYZ_MOM_11) + dx*dx

      xyzMoms(XYZ_MOM_12) = xyzMoms(XYZ_MOM_12) + dx*dy
      xyzMoms(XYZ_MOM_22) = xyzMoms(XYZ_MOM_22) + dy*dy

      xyzMoms(XYZ_MOM_13) = xyzMoms(XYZ_MOM_13) + dx*dz
      xyzMoms(XYZ_MOM_23) = xyzMoms(XYZ_MOM_23) + dy*dz
      xyzMoms(XYZ_MOM_33) = xyzMoms(XYZ_MOM_33) + dz*dz 

      xyzMoms(XYZ_MOM_14) = xyzMoms(XYZ_MOM_14) + wt*dx
      xyzMoms(XYZ_MOM_24) = xyzMoms(XYZ_MOM_24) + wt*dy
      xyzMoms(XYZ_MOM_34) = xyzMoms(XYZ_MOM_34) + wt*dz
      xyzMoms(XYZ_MOM_44) = xyzMoms(XYZ_MOM_44) + wt*wt  
    END DO ! isl

    r11  = SQRT(xyzMoms(XYZ_MOM_11)) 
    ir11 = 1.0_RFREAL/r11

    r12  = ir11*xyzMoms(XYZ_MOM_12)
    r22  = SQRT(xyzMoms(XYZ_MOM_22) - r12*r12)
    ir22 = 1.0_RFREAL/r22

    r13  = ir11*xyzMoms(XYZ_MOM_13)
    r23  = ir22*(xyzMoms(XYZ_MOM_23) -  r12*r13           ) 
    r33  = SQRT(xyzMoms(XYZ_MOM_33)  - (r13*r13 + r23*r23))
    ir33 = 1.0_RFREAL/r33

    r14  = ir11*xyzMoms(XYZ_MOM_14)
    r24  = ir22*(xyzMoms(XYZ_MOM_24) -  r12*r14                     )
    r34  = ir33*(xyzMoms(XYZ_MOM_34) - (r13*r14 + r23*r24          ))
    r44  = SQRT(xyzMoms(XYZ_MOM_44)  - (r14*r14 + r24*r24 + r34*r34))

! ==============================================================================
!   Store weights
! ==============================================================================

    xyzMoms(XYZ_MOM_11) = r11        

    xyzMoms(XYZ_MOM_12) = r12
    xyzMoms(XYZ_MOM_22) = r22

    xyzMoms(XYZ_MOM_13) = r13
    xyzMoms(XYZ_MOM_23) = r23
    xyzMoms(XYZ_MOM_33) = r33 

    xyzMoms(XYZ_MOM_14) = r14
    xyzMoms(XYZ_MOM_24) = r24
    xyzMoms(XYZ_MOM_34) = r34
    xyzMoms(XYZ_MOM_44) = r44
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)    
  
  END SUBROUTINE RFLU_ComputeStencilMoments3D1





! *******************************************************************************
!
! Purpose: Destroy weights for boundary face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch		Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_DestroyWtsBF2C(pRegion,pPatch)

    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    TYPE(t_patch), POINTER :: pPatch      

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid    

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyBWtsF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying boundary face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      DO ifl = 1,pPatch%nBFaces
        DEALLOCATE(pPatch%bf2cs(ifl)%xyzMoms,STAT=errorFlag)  
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2cs%xyzMoms')
        END IF ! global%error
      END DO ! ifl
    END IF ! pPatch%bcType

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsBF2C(pRegion,pPatch)   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying boundary face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsBF2C






! *******************************************************************************
!
! Purpose: Wrapper routine for destroying boundary-face-to-cell weights.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyWtsBF2CWrapper(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyWtsBF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensBFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_DestroyWtsBF2C(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsBF2CWrapper






! *******************************************************************************
!
! Purpose: Destroy weights for cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_DestroyWtsC2C(pRegion)

    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_global), POINTER :: global     
    TYPE(t_grid), POINTER :: pGrid         

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyWtsC2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying cell-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot 
      DEALLOCATE(pGrid%c2cs(icg)%xyzMoms,STAT=errorFlag)                   
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%c2cs%xyzMoms')
      END IF ! global%error        
    END DO ! icg

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsC2C(pRegion)   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying cell-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsC2C









! *******************************************************************************
!
! Purpose: Wrapper routine for destroying cell-to-cell weights.
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

  SUBROUTINE RFLU_DestroyWtsC2CWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyWtsC2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensCells )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_DestroyWtsC2C(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsC2CWrapper








! *******************************************************************************
!
! Purpose: Destroy weights for face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_DestroyWtsF2C(pRegion)

    IMPLICIT NONE
    
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifg
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid    

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyWtsF2C',&
  'RFLU_ModWeights.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying face-to-cell weights...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces
      DEALLOCATE(pGrid%f2cs(ifg)%xyzMoms,STAT=errorFlag)                   
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cs%xyzMoms')
      END IF ! global%error        
    END DO ! ifg

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyWtsF2C(pRegion)   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying face-to-cell weights done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsF2C






! *******************************************************************************
!
! Purpose: Wrapper routine for destroying face-to-cell weights.
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

  SUBROUTINE RFLU_DestroyWtsF2CWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global  
    TYPE(t_mixt_input), POINTER :: pMixtInput   

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyWtsF2CWrapper',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pMixtInput => pRegion%mixtInput 

! ******************************************************************************
!   Call routines to compute weights
! ******************************************************************************

    SELECT CASE ( pMixtInput%stencilDimensFaces )
      CASE ( 1 )               
      CASE ( 2,3 ) 
        CALL RFLU_DestroyWtsF2C(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyWtsF2CWrapper






! *******************************************************************************
!
! Purpose: Nullify weights for boundary face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_NullifyWtsBF2C(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch   
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ifl
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid      

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyWtsBF2C',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      DO ifl = 1,pPatch%nBFaces
        NULLIFY(pPatch%bf2cs(ifl)%xyzMoms)  
      END DO ! ifl
    END IF ! pPatch%bcType

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyWtsBF2C  


  
  
  
  
! *******************************************************************************
!
! Purpose: Nullify weights for cell-to-cell stencil.
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

  SUBROUTINE RFLU_NullifyWtsC2C(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyWtsC2C',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot 
      NULLIFY(pGrid%c2cs(icg)%xyzMoms)                   
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyWtsC2C  
  




! *******************************************************************************
!
! Purpose: Nullify weights for face-to-cell stencil.
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
  
  SUBROUTINE RFLU_NullifyWtsF2C(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ifg
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid          

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyWtsF2C',&
  'RFLU_ModWeights.F90')

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces 
      NULLIFY(pGrid%f2cs(ifg)%xyzMoms)                   
    END DO ! ifg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyWtsF2C  

  
  
  



! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModWeights


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModWeights.F90,v $
! Revision 1.16  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2007/02/27 13:08:27  haselbac
! Enabled 1d computations
!
! Revision 1.13  2006/04/07 16:03:33  haselbac
! Changed computation of bf2c wts to be done patch-wise
!
! Revision 1.12  2006/04/07 15:19:21  haselbac
! Removed tabs
!
! Revision 1.11  2006/04/07 14:52:40  haselbac
! Adapted to changes in stencilDimens params, bug fixes in 1D routine
!
! Revision 1.10  2006/03/09 20:50:56  haselbac
! Bug fix: Removed icg from INTEGER, INTENT(IN)
!
! Revision 1.9  2006/03/09 14:09:22  haselbac
! Wrapperified module bcos of 1D routines
!
! Revision 1.8  2006/01/06 22:14:45  haselbac
! Renamed routines bcos added 1D routines
!
! Revision 1.7  2005/10/05 14:16:51  haselbac
! Added routines for bfaces, clean-up
!
! Revision 1.6  2005/03/09 15:06:03  haselbac
! Added 2d option
!
! Revision 1.5  2004/07/06 15:14:48  haselbac
! Cosmetics only
!                                                
! Revision 1.4  2004/03/18 03:32:39  haselbac                                  
! Bug fix for i4 terms, joined init and computation                            
!
! Revision 1.3  2004/01/22 16:03:59  haselbac                                  
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC   and titan
!
! Revision 1.2  2004/01/13 16:22:32  haselbac                                  
! Changed to inverse-distance weighting                                        
!
! Revision 1.1  2003/12/04 03:28:43  haselbac                                  
! Initial revision                                                             
!
! ******************************************************************************






























