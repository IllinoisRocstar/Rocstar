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
! Purpose: Collection of routines to compute force and moment coefficients.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModForcesMoments.F90,v 1.11 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModForcesMoments

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModForcesMoments.F90,v $ $Revision: 1.11 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ComputeGlobalForcesMoments, &  
            RFLU_ComputeLocalForcesMoments, & 
            RFLU_CreateForcesMoments, & 
            RFLU_CreateGlobalThrustFlags, & 
            RFLU_DestroyForcesMoments, &            
            RFLU_DestroyGlobalThrustFlags, &            
            RFLU_PrintGlobalForcesMoments, & 
            RFLU_SetGlobalThrustFlags, &
            RFLU_WriteGlobalForcesMoments
            
! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_InitForcesMoments, & 
             RFLU_InitGlobalThrustFlags, & 
             RFLU_NullifyForcesMoments, &
             RFLU_NullifyGlobalThrustFlags

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  





! *******************************************************************************
!
! Purpose: Compute global force, moment and mass coefficients.
!
! Description: None.
!
! Input:
!   regions             Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeGlobalForcesMoments(regions)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iComp,iPatch,iPatchGlobal,iReg,iXyz,nVals2,nVals3
    REAL(RFREAL) :: fact
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: globalVals2,localVals2
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: globalVals3,localVals3
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pCoeff     
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => regions(1)%global
    
    CALL RegisterFunction(global,'RFLU_ComputeGlobalForcesMoments',&
  'RFLU_ModForcesMoments.F90')

    nVals2 = (MASS_OUT-MASS_IN+1)*global%nPatches
    nVals3 = (ZCOORD-XCOORD+1)*(COMP_VISC-COMP_MOM+1)*global%nPatches
                
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************
        
    ALLOCATE(globalVals3(XCOORD:ZCOORD,COMP_MOM:COMP_VISC,global%nPatches), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals3')
    END IF ! global%error        
        
    ALLOCATE(localVals3(XCOORD:ZCOORD,COMP_MOM:COMP_VISC,global%nPatches), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals3')
    END IF ! global%error

! ******************************************************************************
!   Compute global force coefficients
! ******************************************************************************

    DO iPatch = 1,global%nPatches
      globalVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
                  
      localVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
    END DO ! iPatch

! ==============================================================================
!   Set local force coefficients 
! ==============================================================================
        
    DO iReg = 1,global%nRegionsLocal   
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        iPatchGlobal = pPatch%iPatchGlobal
             
          localVals3(XCOORD,COMP_MOM ,iPatchGlobal) &
        = localVals3(XCOORD,COMP_MOM ,iPatchGlobal) & 
        + pPatch%forceCoeffs(XCOORD,COMP_MOM ) 

          localVals3(XCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceCoeffs(XCOORD,COMP_PRES)      

          localVals3(XCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceCoeffs(XCOORD,COMP_VISC)      

          localVals3(YCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%forceCoeffs(YCOORD,COMP_MOM )

          localVals3(YCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceCoeffs(YCOORD,COMP_PRES)

          localVals3(YCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceCoeffs(YCOORD,COMP_VISC)

          localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%forceCoeffs(ZCOORD,COMP_MOM )

          localVals3(ZCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceCoeffs(ZCOORD,COMP_PRES)

          localVals3(ZCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceCoeffs(ZCOORD,COMP_VISC)
      END DO ! iPatch
    END DO ! iReg

! ==============================================================================
!   Compute global force coefficients 
! ==============================================================================

    CALL MPI_AllReduce(localVals3,globalVals3,nVals3,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%forceCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = globalVals3(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = globalVals3(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = globalVals3(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = globalVals3(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = globalVals3(YCOORD,COMP_PRES,iPatch)          
        pCoeff(YCOORD,COMP_VISC,iPatch) = globalVals3(YCOORD,COMP_VISC,iPatch)          
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = globalVals3(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = globalVals3(ZCOORD,COMP_PRES,iPatch)          
        pCoeff(ZCOORD,COMP_VISC,iPatch) = globalVals3(ZCOORD,COMP_VISC,iPatch)          
      END DO ! iPatch
    END DO ! iReg
    
! ==============================================================================
!   Normalize force coefficients
! ==============================================================================

    fact = 1.0_RFREAL/global%forceRefArea

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%forceCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = fact*pCoeff(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = fact*pCoeff(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = fact*pCoeff(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = fact*pCoeff(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = fact*pCoeff(YCOORD,COMP_PRES,iPatch)         
        pCoeff(YCOORD,COMP_VISC,iPatch) = fact*pCoeff(YCOORD,COMP_VISC,iPatch)         
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = fact*pCoeff(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = fact*pCoeff(ZCOORD,COMP_PRES,iPatch)       
        pCoeff(ZCOORD,COMP_VISC,iPatch) = fact*pCoeff(ZCOORD,COMP_VISC,iPatch)       
      END DO ! iPatch
    END DO ! iReg 

! ******************************************************************************
!   Compute global vacuum force coefficients
! ******************************************************************************

    DO iPatch = 1,global%nPatches
      globalVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
                  
      localVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
    END DO ! iPatch

! ==============================================================================
!   Set local vacuum force coefficients 
! ==============================================================================
        
    DO iReg = 1,global%nRegionsLocal   
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        iPatchGlobal = pPatch%iPatchGlobal
             
          localVals3(XCOORD,COMP_MOM ,iPatchGlobal) &
        = localVals3(XCOORD,COMP_MOM ,iPatchGlobal) & 
        + pPatch%forceVacCoeffs(XCOORD,COMP_MOM ) 

          localVals3(XCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceVacCoeffs(XCOORD,COMP_PRES)      

          localVals3(XCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceVacCoeffs(XCOORD,COMP_VISC)      

          localVals3(YCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%forceVacCoeffs(YCOORD,COMP_MOM )

          localVals3(YCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceVacCoeffs(YCOORD,COMP_PRES)

          localVals3(YCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceVacCoeffs(YCOORD,COMP_VISC)

          localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%forceVacCoeffs(ZCOORD,COMP_MOM )

          localVals3(ZCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%forceVacCoeffs(ZCOORD,COMP_PRES)

          localVals3(ZCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%forceVacCoeffs(ZCOORD,COMP_VISC)
      END DO ! iPatch
    END DO ! iReg

! ==============================================================================
!   Compute global vacuum force coefficients 
! ==============================================================================

    CALL MPI_AllReduce(localVals3,globalVals3,nVals3,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%forceVacCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = globalVals3(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = globalVals3(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = globalVals3(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = globalVals3(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = globalVals3(YCOORD,COMP_PRES,iPatch)          
        pCoeff(YCOORD,COMP_VISC,iPatch) = globalVals3(YCOORD,COMP_VISC,iPatch)          
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = globalVals3(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = globalVals3(ZCOORD,COMP_PRES,iPatch)          
        pCoeff(ZCOORD,COMP_VISC,iPatch) = globalVals3(ZCOORD,COMP_VISC,iPatch)          
      END DO ! iPatch
    END DO ! iReg
    
! ==============================================================================
!   Normalize vacuum force coefficients
! ==============================================================================

    fact = 1.0_RFREAL/global%forceRefArea

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%forceVacCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = fact*pCoeff(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = fact*pCoeff(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = fact*pCoeff(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = fact*pCoeff(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = fact*pCoeff(YCOORD,COMP_PRES,iPatch)         
        pCoeff(YCOORD,COMP_VISC,iPatch) = fact*pCoeff(YCOORD,COMP_VISC,iPatch)         
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = fact*pCoeff(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = fact*pCoeff(ZCOORD,COMP_PRES,iPatch)       
        pCoeff(ZCOORD,COMP_VISC,iPatch) = fact*pCoeff(ZCOORD,COMP_VISC,iPatch)       
      END DO ! iPatch
    END DO ! iReg 

! ******************************************************************************
!   Compute global moment coefficients
! ******************************************************************************

    DO iPatch = 1,global%nPatches
      globalVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      globalVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      globalVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
                  
      localVals3(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      localVals3(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL      
      localVals3(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL      
    END DO ! iPatch

! ==============================================================================
!   Set local moment coefficients
! ==============================================================================
        
    DO iReg = 1,global%nRegionsLocal   
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        iPatchGlobal = pPatch%iPatchGlobal
             
          localVals3(XCOORD,COMP_MOM ,iPatchGlobal) &
        = localVals3(XCOORD,COMP_MOM ,iPatchGlobal) & 
        + pPatch%momentCoeffs(XCOORD,COMP_MOM ) 

          localVals3(XCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%momentCoeffs(XCOORD,COMP_PRES)      

          localVals3(XCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(XCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%momentCoeffs(XCOORD,COMP_VISC)      

          localVals3(YCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%momentCoeffs(YCOORD,COMP_MOM )

          localVals3(YCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%momentCoeffs(YCOORD,COMP_PRES)

          localVals3(YCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(YCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%momentCoeffs(YCOORD,COMP_VISC)

          localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_MOM ,iPatchGlobal) &
        + pPatch%momentCoeffs(ZCOORD,COMP_MOM )

          localVals3(ZCOORD,COMP_PRES,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_PRES,iPatchGlobal) &
        + pPatch%momentCoeffs(ZCOORD,COMP_PRES)

          localVals3(ZCOORD,COMP_VISC,iPatchGlobal) & 
        = localVals3(ZCOORD,COMP_VISC,iPatchGlobal) &
        + pPatch%momentCoeffs(ZCOORD,COMP_VISC)
      END DO ! iPatch
    END DO ! iReg

! ==============================================================================
!   Compute global moment coefficients
! ==============================================================================

    CALL MPI_AllReduce(localVals3,globalVals3,nVals3,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%momentCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = globalVals3(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = globalVals3(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = globalVals3(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = globalVals3(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = globalVals3(YCOORD,COMP_PRES,iPatch)          
        pCoeff(YCOORD,COMP_VISC,iPatch) = globalVals3(YCOORD,COMP_VISC,iPatch)          
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = globalVals3(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = globalVals3(ZCOORD,COMP_PRES,iPatch)          
        pCoeff(ZCOORD,COMP_VISC,iPatch) = globalVals3(ZCOORD,COMP_VISC,iPatch)          
      END DO ! iPatch
    END DO ! iReg
    
! ==============================================================================
!   Normalize moment coefficients
! ==============================================================================

    fact = 1.0_RFREAL/(global%forceRefArea*global%forceRefLength)

    DO iReg = 1,global%nRegionsLocal   
      pCoeff => regions(iReg)%momentCoeffsGlobal
        
      DO iPatch = 1,pRegion%global%nPatches
        pCoeff(XCOORD,COMP_MOM ,iPatch) = fact*pCoeff(XCOORD,COMP_MOM ,iPatch)
        pCoeff(XCOORD,COMP_PRES,iPatch) = fact*pCoeff(XCOORD,COMP_PRES,iPatch)
        pCoeff(XCOORD,COMP_VISC,iPatch) = fact*pCoeff(XCOORD,COMP_VISC,iPatch)
        pCoeff(YCOORD,COMP_MOM ,iPatch) = fact*pCoeff(YCOORD,COMP_MOM ,iPatch)
        pCoeff(YCOORD,COMP_PRES,iPatch) = fact*pCoeff(YCOORD,COMP_PRES,iPatch)          
        pCoeff(YCOORD,COMP_VISC,iPatch) = fact*pCoeff(YCOORD,COMP_VISC,iPatch)          
        pCoeff(ZCOORD,COMP_MOM ,iPatch) = fact*pCoeff(ZCOORD,COMP_MOM ,iPatch)
        pCoeff(ZCOORD,COMP_PRES,iPatch) = fact*pCoeff(ZCOORD,COMP_PRES,iPatch)         
        pCoeff(ZCOORD,COMP_VISC,iPatch) = fact*pCoeff(ZCOORD,COMP_VISC,iPatch)         
      END DO ! iPatch
    END DO ! iReg

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************
        
    DEALLOCATE(globalVals3,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals3')
    END IF ! global%error        
        
    DEALLOCATE(localVals3,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals3')
    END IF ! global%error

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(globalVals2(MASS_IN:MASS_OUT,global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals2')
    END IF ! global%error        
        
    ALLOCATE(localVals2(MASS_IN:MASS_OUT,global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals2')
    END IF ! global%error

! ******************************************************************************
!   Compute global mass coefficients
! ******************************************************************************

    DO iPatch = 1,global%nPatches
      globalVals2(MASS_IN ,iPatch) = 0.0_RFREAL
      globalVals2(MASS_OUT,iPatch) = 0.0_RFREAL

      localVals2(MASS_IN ,iPatch)  = 0.0_RFREAL
      localVals2(MASS_OUT,iPatch)  = 0.0_RFREAL
    END DO ! iPatch

! ==============================================================================
!   Set local mass coefficients 
! ==============================================================================
        
    DO iReg = 1,global%nRegionsLocal   
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        iPatchGlobal = pPatch%iPatchGlobal
             
        localVals2(MASS_IN ,iPatchGlobal) = localVals2(MASS_IN ,iPatchGlobal) &
                                          + pPatch%massCoeffs(MASS_IN)
        localVals2(MASS_OUT,iPatchGlobal) = localVals2(MASS_OUT,iPatchGlobal) &
                                          + pPatch%massCoeffs(MASS_OUT)
      END DO ! iPatch
    END DO ! iReg

! ==============================================================================
!   Compute global mass coefficients 
! ==============================================================================

    CALL MPI_AllReduce(localVals2,globalVals2,nVals2,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal   
      DO iPatch = 1,pRegion%global%nPatches
        regions(iReg)%massCoeffsGlobal(MASS_IN,iPatch)  &
          = globalVals2(MASS_IN,iPatch)
        regions(iReg)%massCoeffsGlobal(MASS_OUT,iPatch) &
          = globalVals2(MASS_OUT,iPatch)
      END DO ! iPatch
    END DO ! iReg
    
! ==============================================================================
!   Normalize mass coefficients
! ==============================================================================

    fact = 1.0_RFREAL/global%forceRefArea

    DO iReg = 1,global%nRegionsLocal   
      DO iPatch = 1,pRegion%global%nPatches
        regions(iReg)%massCoeffsGlobal(MASS_IN,iPatch)  &
          = fact*regions(iReg)%massCoeffsGlobal(MASS_IN,iPatch)
        regions(iReg)%massCoeffsGlobal(MASS_OUT,iPatch) &
          = fact*regions(iReg)%massCoeffsGlobal(MASS_OUT,iPatch)
      END DO ! iPatch
    END DO ! iReg

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************
     
    DEALLOCATE(globalVals2,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals2')
    END IF ! global%error        
        
    DEALLOCATE(localVals2,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals2')
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGlobalForcesMoments








! *******************************************************************************
!
! Purpose: Compute local force, moment and coefficients.
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

  SUBROUTINE RFLU_ComputeLocalForcesMoments(pRegion)

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
    
    INTEGER :: ifl,iPatch
    REAL(RFREAL) :: cfx,cfy,cfz,ch,cmass,cmomx,cmomy,cmomz,cp,cpref,fmx,fpx, &
                    fpx_vac,fmy,fpy,fpy_vac,fmz,fpz,fpz_vac,fvx,fvy,fvz,mc_in, &
                    mc_out,mpx,mpy,mpz,mvx,mvy,mvz,nm,nx,ny,nz,xc,xRef,yc, &
                    yRef,zc,zRef
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ComputeLocalForcesMoments',&
  'RFLU_ModForcesMoments.F90')
    
    pGrid => pRegion%grid
            
    xRef = global%forceRefXCoord
    yRef = global%forceRefYCoord
    zRef = global%forceRefZCoord

    cpref = global%refPressure/ &
            (0.5_RFREAL*global%refDensity*global%refVelocity*global%refVelocity)

! ******************************************************************************
!   Loop over patches
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
            
      fpx = 0.0_RFREAL
      fpy = 0.0_RFREAL
      fpz = 0.0_RFREAL
                        
      fpx_vac = 0.0_RFREAL
      fpy_vac = 0.0_RFREAL
      fpz_vac = 0.0_RFREAL
                        
      fvx = 0.0_RFREAL
      fvy = 0.0_RFREAL
      fvz = 0.0_RFREAL
      
      fmx = 0.0_RFREAL
      fmy = 0.0_RFREAL
      fmz = 0.0_RFREAL
      
      mpx = 0.0_RFREAL
      mpy = 0.0_RFREAL
      mpz = 0.0_RFREAL
                     
      mvx = 0.0_RFREAL
      mvy = 0.0_RFREAL
      mvz = 0.0_RFREAL      

      mc_in  = 0.0_RFREAL                   
      mc_out = 0.0_RFREAL                   
 
! ==============================================================================
!     Loop over faces 
! ==============================================================================                    
                                    
      DO ifl = 1,pPatch%nBFaces
      
! ------------------------------------------------------------------------------
!       Get geometry
! ------------------------------------------------------------------------------      
      
        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)                
        nm = pPatch%fn(XYZMAG,ifl)
      
        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)                
      
! ------------------------------------------------------------------------------
!       Get coefficients
! ------------------------------------------------------------------------------      
      
        cp  = pPatch%cp(ifl)
        cfx = pPatch%cf(XCOORD,ifl)
        cfy = pPatch%cf(YCOORD,ifl)
        cfz = pPatch%cf(ZCOORD,ifl)
        ch  = pPatch%ch(ifl)                

        cmass = pPatch%cmass(ifl)                
        cmomx = pPatch%cmom(XCOORD,ifl)                
        cmomy = pPatch%cmom(YCOORD,ifl)                
        cmomz = pPatch%cmom(ZCOORD,ifl)                

! ------------------------------------------------------------------------------
!       Compute contributions to force and moment coefficients
! ------------------------------------------------------------------------------      
      
        fpx = fpx + cp*nx*nm
        fpy = fpy + cp*ny*nm
        fpz = fpz + cp*nz*nm

        fpx_vac = fpx_vac + (cp+cpref)*nx*nm
        fpy_vac = fpy_vac + (cp+cpref)*ny*nm
        fpz_vac = fpz_vac + (cp+cpref)*nz*nm

        fvx = fvx + cfx*nm
        fvy = fvy + cfy*nm
        fvz = fvz + cfz*nm
                        
        fmx = fmx + cmomx*nm
        fmy = fmy + cmomy*nm
        fmz = fmz + cmomz*nm
                        
        mpx = mpx - cp*(ny*(zc - zRef) + nz*(yc - yRef))*nm
        mpy = mpy + cp*(nx*(zc - zRef) - nz*(xc - xRef))*nm
        mpz = mpz + cp*(ny*(xc - xRef) - nx*(yc - yRef))*nm                 

        mvx = mvx - (cfy*(zc - zRef) + cfz*(yc - yRef))*nm
        mvy = mvy + (cfx*(zc - zRef) - cfz*(xc - xRef))*nm
        mvz = mvz + (cfy*(xc - xRef) - cfx*(yc - yRef))*nm                 

        mc_in  = mc_in  - 0.5_RFREAL*(cmass-ABS(cmass))*nm
        mc_out = mc_out + 0.5_RFREAL*(cmass+ABS(cmass))*nm
      END DO ! ifl
          
! ==============================================================================
!     Normalize and store coefficients
! ==============================================================================        
                    
      pPatch%forceCoeffs(XCOORD,COMP_MOM ) = fmx
      pPatch%forceCoeffs(YCOORD,COMP_MOM ) = fmy
      pPatch%forceCoeffs(ZCOORD,COMP_MOM ) = fmz
       
      pPatch%forceCoeffs(XCOORD,COMP_PRES) = fpx
      pPatch%forceCoeffs(YCOORD,COMP_PRES) = fpy
      pPatch%forceCoeffs(ZCOORD,COMP_PRES) = fpz   

      pPatch%forceCoeffs(XCOORD,COMP_VISC) = fvx
      pPatch%forceCoeffs(YCOORD,COMP_VISC) = fvy
      pPatch%forceCoeffs(ZCOORD,COMP_VISC) = fvz   

      pPatch%forceVacCoeffs(XCOORD,COMP_MOM ) = fmx 
      pPatch%forceVacCoeffs(YCOORD,COMP_MOM ) = fmy 
      pPatch%forceVacCoeffs(ZCOORD,COMP_MOM ) = fmz
       
      pPatch%forceVacCoeffs(XCOORD,COMP_PRES) = fpx_vac
      pPatch%forceVacCoeffs(YCOORD,COMP_PRES) = fpx_vac
      pPatch%forceVacCoeffs(ZCOORD,COMP_PRES) = fpx_vac   

      pPatch%forceVacCoeffs(XCOORD,COMP_VISC) = fvx
      pPatch%forceVacCoeffs(YCOORD,COMP_VISC) = fvy
      pPatch%forceVacCoeffs(ZCOORD,COMP_VISC) = fvz   

      pPatch%momentCoeffs(XCOORD,COMP_PRES) = mpx
      pPatch%momentCoeffs(YCOORD,COMP_PRES) = mpy
      pPatch%momentCoeffs(ZCOORD,COMP_PRES) = mpz
       
      pPatch%momentCoeffs(XCOORD,COMP_VISC) = mvx
      pPatch%momentCoeffs(YCOORD,COMP_VISC) = mvy
      pPatch%momentCoeffs(ZCOORD,COMP_VISC) = mvz   
       
      pPatch%massCoeffs(MASS_IN)  = mc_in
      pPatch%massCoeffs(MASS_OUT) = mc_out
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeLocalForcesMoments





  


! *******************************************************************************
!
! Purpose: Create force,moment and mass coefficients.
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

  SUBROUTINE RFLU_CreateForcesMoments(pRegion)

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
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_CreateForcesMoments',&
  'RFLU_ModForcesMoments.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      ALLOCATE(pPatch%forceCoeffs(XCOORD:ZCOORD,COMP_MOM :COMP_VISC), & 
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%forceCoeffs')
      END IF ! global%error  
            
      ALLOCATE(pPatch%forceVacCoeffs(XCOORD:ZCOORD,COMP_MOM :COMP_VISC), & 
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%forceVacCoeffs')
      END IF ! global%error  
            
      ALLOCATE(pPatch%momentCoeffs(XCOORD:ZCOORD,COMP_MOM :COMP_VISC), & 
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%momentCoeffs')
      END IF ! global%error                    
      
      ALLOCATE(pPatch%massCoeffs(MASS_IN:MASS_OUT),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%massCoeffs')
      END IF ! global%error  
    END DO ! iPatch

    ALLOCATE(pRegion%forceCoeffsGlobal(XCOORD:ZCOORD, &
             COMP_MOM:COMP_VISC,global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%forceCoeffsGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%forceVacCoeffsGlobal(XCOORD:ZCOORD, &
             COMP_MOM:COMP_VISC,global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%forceVacCoeffsGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%momentCoeffsGlobal(XCOORD:ZCOORD, &
             COMP_MOM:COMP_VISC,global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                     'pRegion%momentCoeffsGlobal')
    END IF ! global%error

    ALLOCATE(pRegion%massCoeffsGlobal(MASS_IN:MASS_OUT, &
             global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%massCoeffsGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%specImpulseGlobal(XCOORD:ZCOORD, &
             global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%specImpulseGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%specImpulseVacGlobal(XCOORD:ZCOORD, &
             global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%specImpulseVacGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%thrustGlobal(XCOORD:ZCOORD, &
             global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%thrustGlobal')
    END IF ! global%error  
          
    ALLOCATE(pRegion%thrustVacGlobal(XCOORD:ZCOORD, &
             global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%thrustVacGlobal')
    END IF ! global%error  
          
! ******************************************************************************
!   Initialize memory
! ******************************************************************************

    CALL RFLU_InitForcesMoments(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateForcesMoments
 








! *******************************************************************************
!
! Purpose: Create global thrust flags.
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

  SUBROUTINE RFLU_CreateGlobalThrustFlags(pRegion)

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
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_CreateGlobalThrustFlags',&
  'RFLU_ModForcesMoments.F90')
    
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    ALLOCATE(pRegion%thrustFlagsGlobal(global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%thrustFlagsGlobal')
    END IF ! global%error  
          
! ******************************************************************************
!   Initialize memory
! ******************************************************************************

    CALL RFLU_InitGlobalThrustFlags(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateGlobalThrustFlags








! *******************************************************************************
!
! Purpose: Destroy force, moment and mass coefficients.
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

  SUBROUTINE RFLU_DestroyForcesMoments(pRegion)

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
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_DestroyForcesMoments',&
  'RFLU_ModForcesMoments.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%forceCoeffs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%forceCoeffs')
      END IF ! global%error  
            
      DEALLOCATE(pPatch%forceVacCoeffs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%forceVacCoeffs')
      END IF ! global%error  
            
      DEALLOCATE(pPatch%momentCoeffs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%momentCoeffs')
      END IF ! global%error                    
            
      DEALLOCATE(pPatch%massCoeffs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%massCoeffs')
      END IF ! global%error  
    END DO ! iPatch

    DEALLOCATE(pRegion%forceCoeffsGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                     'pRegion%forceCoeffsGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%forceVacCoeffsGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                     'pRegion%forceVacCoeffsGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%momentCoeffsGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                     'pRegion%momentCoeffsGlobal')
    END IF ! global%error            

    DEALLOCATE(pRegion%massCoeffsGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%massCoeffsGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%specImpulseGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%specImpulseGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%specImpulseVacGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%specImpulseVacGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%thrustGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%thrustGlobal')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%thrustVacGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%thrustVacGlobal')
    END IF ! global%error  
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyForcesMoments(pRegion)
    
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyForcesMoments









! *******************************************************************************
!
! Purpose: Destroy global thrust flags 
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

  SUBROUTINE RFLU_DestroyGlobalThrustFlags(pRegion)

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
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_DestroyGlobalThrustFlag',&
  'RFLU_ModForcesMoments.F90')
    
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DEALLOCATE(pRegion%thrustFlagsGlobal,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%thrustFlagsGlobal')
    END IF ! global%error  
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************
        
    CALL RFLU_NullifyGlobalThrustFlags(pRegion)
          
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyGlobalThrustFlags







! *******************************************************************************
!
! Purpose: Initialize force, moment and mass coefficients.
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

  SUBROUTINE RFLU_InitForcesMoments(pRegion)

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
    
    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_InitForcesMoments',&
  'RFLU_ModForcesMoments.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      pPatch%forceCoeffs(XCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceCoeffs(XCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%forceCoeffs(XCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%forceCoeffs(YCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceCoeffs(YCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%forceCoeffs(YCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%forceCoeffs(ZCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceCoeffs(ZCOORD,COMP_PRES) = 0.0_RFREAL                  
      pPatch%forceCoeffs(ZCOORD,COMP_VISC) = 0.0_RFREAL                  

      pPatch%forceVacCoeffs(XCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceVacCoeffs(XCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%forceVacCoeffs(XCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%forceVacCoeffs(YCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceVacCoeffs(YCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%forceVacCoeffs(YCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%forceVacCoeffs(ZCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%forceVacCoeffs(ZCOORD,COMP_PRES) = 0.0_RFREAL                  
      pPatch%forceVacCoeffs(ZCOORD,COMP_VISC) = 0.0_RFREAL                  

      pPatch%momentCoeffs(XCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%momentCoeffs(XCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%momentCoeffs(XCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%momentCoeffs(YCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%momentCoeffs(YCOORD,COMP_PRES) = 0.0_RFREAL      
      pPatch%momentCoeffs(YCOORD,COMP_VISC) = 0.0_RFREAL      
      pPatch%momentCoeffs(ZCOORD,COMP_MOM ) = 0.0_RFREAL
      pPatch%momentCoeffs(ZCOORD,COMP_PRES) = 0.0_RFREAL                  
      pPatch%momentCoeffs(ZCOORD,COMP_VISC) = 0.0_RFREAL                  

      pPatch%massCoeffs(MASS_IN)  = 0.0_RFREAL
      pPatch%massCoeffs(MASS_OUT) = 0.0_RFREAL                  
    END DO ! iPatch

    DO iPatch = 1,pRegion%global%nPatches
      pRegion%forceCoeffsGlobal(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL
      pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL
      pRegion%forceCoeffsGlobal(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL 
      pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL 
      pRegion%forceCoeffsGlobal(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL 
      pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL 
      pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL 
      
      pRegion%forceVacCoeffsGlobal(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%forceVacCoeffsGlobal(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL
      pRegion%forceVacCoeffsGlobal(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL
      pRegion%forceVacCoeffsGlobal(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%forceVacCoeffsGlobal(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL 
      pRegion%forceVacCoeffsGlobal(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL 
      pRegion%forceVacCoeffsGlobal(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL 
      pRegion%forceVacCoeffsGlobal(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL 
      pRegion%forceVacCoeffsGlobal(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL 
      
      pRegion%momentCoeffsGlobal(XCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch) = 0.0_RFREAL
      pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch) = 0.0_RFREAL
      pRegion%momentCoeffsGlobal(YCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL
      pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch) = 0.0_RFREAL 
      pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch) = 0.0_RFREAL 
      pRegion%momentCoeffsGlobal(ZCOORD,COMP_MOM ,iPatch) = 0.0_RFREAL 
      pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch) = 0.0_RFREAL       
      pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch) = 0.0_RFREAL       
      
      pRegion%massCoeffsGlobal(MASS_IN ,iPatch) = 0.0_RFREAL
      pRegion%massCoeffsGlobal(MASS_OUT,iPatch) = 0.0_RFREAL
      
      pRegion%specImpulseGlobal(XCOORD,iPatch) = 0.0_RFREAL
      pRegion%specImpulseGlobal(YCOORD,iPatch) = 0.0_RFREAL
      pRegion%specImpulseGlobal(ZCOORD,iPatch) = 0.0_RFREAL
      
      pRegion%specImpulseVacGlobal(XCOORD,iPatch)  = 0.0_RFREAL
      pRegion%specImpulseVacGlobal(YCOORD,iPatch)  = 0.0_RFREAL
      pRegion%specImpulseVacGlobal(ZCOORD,iPatch)  = 0.0_RFREAL
      
      pRegion%thrustGlobal(XCOORD,iPatch) = 0.0_RFREAL
      pRegion%thrustGlobal(YCOORD,iPatch) = 0.0_RFREAL
      pRegion%thrustGlobal(ZCOORD,iPatch) = 0.0_RFREAL
      
      pRegion%thrustVacGlobal(XCOORD,iPatch)  = 0.0_RFREAL
      pRegion%thrustVacGlobal(YCOORD,iPatch)  = 0.0_RFREAL
      pRegion%thrustVacGlobal(ZCOORD,iPatch)  = 0.0_RFREAL
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InitForcesMoments







! *******************************************************************************
!
! Purpose: Initialize global thrust flags.
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

  SUBROUTINE RFLU_InitGlobalThrustFlags(pRegion)

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

    INTEGER :: iPatch    
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_InitGlobalThrustFlags',&
  'RFLU_ModForcesMoments.F90')
    
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pRegion%global%nPatches
      pRegion%thrustFlagsGlobal(iPatch) = .FALSE.
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InitGlobalThrustFlags








! *******************************************************************************
!
! Purpose: Nullify force, moment and mass coefficients.
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

  SUBROUTINE RFLU_NullifyForcesMoments(pRegion)

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
    
    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_NullifyForcesMoments',&
  'RFLU_ModForcesMoments.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      NULLIFY(pPatch%forceCoeffs)
      NULLIFY(pPatch%forceVacCoeffs)
      NULLIFY(pPatch%momentCoeffs)
      NULLIFY(pPatch%massCoeffs)
    END DO ! iPatch

    NULLIFY(pRegion%forceCoeffsGlobal)
    NULLIFY(pRegion%forceVacCoeffsGlobal)
    NULLIFY(pRegion%momentCoeffsGlobal)
    NULLIFY(pRegion%massCoeffsGlobal)
    NULLIFY(pRegion%specImpulseGlobal)
    NULLIFY(pRegion%specImpulseVacGlobal)
    NULLIFY(pRegion%thrustGlobal)
    NULLIFY(pRegion%thrustVacGlobal)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyForcesMoments







! *******************************************************************************
!
! Purpose: Nullify global thrust flags 
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

  SUBROUTINE RFLU_NullifyGlobalThrustFlags(pRegion)

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
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_NullifyGlobalThrustFlag',&
  'RFLU_ModForcesMoments.F90')
    
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    NULLIFY(pRegion%thrustFlagsGlobal)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyGlobalThrustFlags







! *******************************************************************************
!
! Purpose: Print global force and moment coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!
! ******************************************************************************

  SUBROUTINE RFLU_PrintGlobalForcesMoments(pRegion)

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
    
    INTEGER :: iPatch 
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_PrintGlobalForcesMoments',&
  'RFLU_ModForcesMoments.F90')
                
! ******************************************************************************
!   Print global force coefficients
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing global force and '// &
                               'moment coefficients...' 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Force coefficients:'                             
                               
      DO iPatch = 1,global%nPatches
        WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Patch:',iPatch
        WRITE(STDOUT,'(A,7X,A,7X,A,6X,A,8X,A)') SOLVER_NAME, & 
              'Component','Pressure','Viscous','Total'
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'x-direction:', & 
              pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
              pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
              pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
              pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch)
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'y-direction:', & 
              pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
              pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
              pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
              pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch)
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'z-direction:', & 
              pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
              pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
              pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
              pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch)                  
      END DO ! iPatch       
      
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Moment coefficients:'                             
                               
      DO iPatch = 1,global%nPatches
        WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Patch:',iPatch
        WRITE(STDOUT,'(A,7X,A,7X,A,6X,A,8X,A)') SOLVER_NAME, & 
              'Component','Pressure','Viscous','Total'
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'x-direction:', & 
              pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
              pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
              pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
              pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch)
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'y-direction:', & 
              pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
              pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
              pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
              pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch)
        WRITE(STDOUT,'(A,7X,A,3(1X,E13.6))') SOLVER_NAME,'z-direction:', & 
              pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
              pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
              pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
              pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch)               
      END DO ! iPatch        
      
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing global force and '// & 
                               'moment coefficients done.'             
    END IF ! global%myProcid

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PrintGlobalForcesMoments







! *******************************************************************************
!
! Purpose: Set global thrust flags.
!
! Description: None.
!
! Input:
!   regions             Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetGlobalThrustFlags(regions)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch,iPatchGlobal,iReg,nVals
    INTEGER, DIMENSION(:), ALLOCATABLE :: globalVals,localVals
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => regions(1)%global
    
    CALL RegisterFunction(global,'RFLU_SetGlobalThrustFlags',&
  'RFLU_ModForcesMoments.F90')

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVals = global%nPatches
                
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************
        
    ALLOCATE(globalVals(global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals')
    END IF ! global%error        
        
    ALLOCATE(localVals(global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals')
    END IF ! global%error

! ******************************************************************************
!   Compute global thrust coefficients flags
! ******************************************************************************

    DO iPatch = 1,global%nPatches
      globalVals(iPatch)  = 0 ! NOTE must be zero
      localVals(iPatch)   = 0
    END DO ! iPatch

! ==============================================================================
!   Set local thrust coefficients flags
! ==============================================================================
        
    DO iReg = 1,global%nRegionsLocal   
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        iPatchGlobal = pPatch%iPatchGlobal
             
        IF ( pPatch%thrustFlag .EQV. .TRUE. ) THEN
          localVals(iPatchGlobal) = localVals(iPatchGlobal) + 1
        END IF ! pPatch%thrustFlag
      END DO ! iPatch
    END DO ! iReg

! ==============================================================================
!   Compute global thrust coefficients flags
! ==============================================================================

    CALL MPI_AllReduce(localVals,globalVals,nVals,MPI_INTEGER,MPI_SUM, &
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal   
      DO iPatch = 1,pRegion%global%nPatches
        IF ( globalVals(iPatch) > 0 ) THEN
          regions(iReg)%thrustFlagsGlobal(iPatch) = .TRUE. 
        END IF ! pPatch%thrustFlag
      END DO ! iPatch
    END DO ! iReg
    
! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************
        
    DEALLOCATE(globalVals,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals')
    END IF ! global%error        
        
    DEALLOCATE(localVals,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals')
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SetGlobalThrustFlags








! *******************************************************************************
!
! Purpose: Write global force and moment coefficients to file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!   2. At present, each patch is written to separate file. This may change in
!      future with definition of superpatches which are treated together. 
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteGlobalForcesMoments(pRegion)

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
    
    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile,iPatch 
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WriteGlobalForcesMoments',&
  'RFLU_ModForcesMoments.F90')

! ******************************************************************************
!   Increment counter
! ******************************************************************************

    global%forceWriteCntr = global%forceWriteCntr + 1
                
! ******************************************************************************
!   Write global force coefficients
! ******************************************************************************

    iFile = IF_FORMOM

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_LOW ) THEN   
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing global force and '// &
                                 'moment coefficients...' 
      END IF ! global%verbLevel

! ==============================================================================
!     Loop over patches
! ==============================================================================

      DO iPatch = 1,global%nPatches
        IF ( global%verbLevel > VERBOSE_LOW ) THEN   
          WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Patch:',iPatch 
        END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!       Build file name
! ------------------------------------------------------------------------------

        WRITE(iFileName,'(A,I4.4)') TRIM(global%outDir)// & 
                                    TRIM(global%casename)//'.fom_',iPatch

! ------------------------------------------------------------------------------
!       Open file
! ------------------------------------------------------------------------------

        IF ( global%restartFromScratch .EQV. .FALSE. ) THEN
          INQUIRE(FILE=iFileName,EXIST=fileExists)

          IF ( fileExists .EQV. .TRUE. ) THEN 
            OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
                 POSITION='APPEND',IOSTAT=errorFlag)
          ELSE 
            OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='NEW', &
                 IOSTAT=errorFlag)
          END IF ! fileExists
        ELSE
          IF ( global%forceWriteCntr == 1 ) THEN 
            OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', &
                 IOSTAT=errorFlag)
          ELSE 
            OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
                 POSITION='APPEND',IOSTAT=errorFlag)          
          END IF ! global%forceWriteCntr
        END IF ! global      

        global%error = errorFlag

        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '// & 
                         TRIM(iFileName))
        END IF ! global%error

! ------------------------------------------------------------------------------
!       Write to file
! ------------------------------------------------------------------------------

        IF ( global%flowType == FLOW_STEADY ) THEN                               
          WRITE(iFile,'(I6,18(1X,E13.6))') global%currentIter, & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch), &

             pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
             pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &

             pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), &

             pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch), &

             pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
             pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch), & 

             pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch) 
        ELSE 
          WRITE(iFile,'(1PE12.5,18(1X,E13.6))') global%currentTime, & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(XCOORD,COMP_VISC,iPatch), &

             pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
             pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &

             pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
             pRegion%forceCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), &

             pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch), & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(XCOORD,COMP_VISC,iPatch), &

             pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch), &
             pRegion%momentCoeffsGlobal(YCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(YCOORD,COMP_VISC,iPatch), & 

             pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch), & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch), & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_PRES,iPatch)+ & 
             pRegion%momentCoeffsGlobal(ZCOORD,COMP_VISC,iPatch)         
         END IF ! global%flowType                                     

! ------------------------------------------------------------------------------
!       Close file
! ------------------------------------------------------------------------------

        CLOSE(iFile,IOSTAT=errorFlag)
        IF (global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                         TRIM(iFileName))
        END IF ! global%error
      END DO ! iPatch        
            
      IF ( global%verbLevel > VERBOSE_LOW ) THEN       
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing global force and '// & 
                                 'moment coefficients done.'             
      END IF ! global%verbLevel
    END IF ! global%myProcid

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteGlobalForcesMoments






END MODULE RFLU_ModForcesMoments

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModForcesMoments.F90,v $
! Revision 1.11  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2007/12/19 19:13:42  rfiedler
! Increased verbosity level required for writing to screen dump to > VERBOSE_LOW.
!
! Revision 1.8  2006/11/02 21:08:34  mparmar
! Bug fix in storing local force-coeffs and normalizing global mass-coeffs
!
! Revision 1.7  2006/10/20 21:18:05  mparmar
! Added code for thrust and specific impulse computation
!
! Revision 1.6  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.5  2005/04/15 15:06:53  haselbac
! Converted to MPI
!
! Revision 1.4  2005/01/05 01:43:32  haselbac
! Added init of global coeffs, now also defined for all regions
!
! Revision 1.3  2004/07/06 15:14:40  haselbac
! Cosmetics only
!                               
! Revision 1.2  2004/06/25 20:07:47  haselbac                       
! Various improvements, fixed bug in file output format statements  
!
! Revision 1.1  2004/06/16 20:00:59  haselbac                       
! Initial revision                                                  
!
! ******************************************************************************



















