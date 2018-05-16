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
! Purpose: Collection of routines for patch operations.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPatchUtils.F90,v 1.7 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPatchUtils

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
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPatchUtils.F90,v $ $Revision: 1.7 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_BuildPatchNeighborMaps, &
            RFLU_CheckPatchBcConsistency, &
            RFLU_ComputePatchNormalsGlobal, & 
            RFLU_ComputePatchNormalsLocal, & 
            RFLU_CreatePatchNeighborMaps, &
            RFLU_DestroyPatchNeighborMaps, & 
            RFLU_GetPatchNormalDirection 

! ==============================================================================
! Private functions
! ==============================================================================



! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
    
  



! ******************************************************************************
!
! Purpose: Build patch-neighbor maps, which indicate whether any two patches 
!   share vertices.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_BuildPatchNeighborMaps(pRegion)

  USE ModSortSearch

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

  INTEGER :: errorFlag,iLoc,iPatch,iPatch2,ivg,ivl
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch,pPatch2

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BuildPatchNeighborMaps',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building patch-neighbor maps...' 
  END IF ! global%myProcid
 
  pGrid => pRegion%grid
 
! ******************************************************************************
! Loop over patches
! ******************************************************************************

  IF ( pGrid%nPatches > 0 ) THEN 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      DO iPatch2 = iPatch+1,pGrid%nPatches
        pPatch2 => pRegion%patches(iPatch2)
      
        ivlLoop: DO ivl = 1,pPatch%nBVert
          ivg = pPatch%bv(ivl)

          CALL BinarySearchInteger(pPatch2%bv(1:pPatch2%nBVert),pPatch2%nBVert, & 
                                   ivg,iLoc)

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            pPatch%nbMap(pPatch2%iPatchGlobal) = .TRUE.

            EXIT ivlLoop
          END IF ! iLoc        
        END DO ivlLoop
      END DO ! iPatch2
    END DO ! iPatch  
    
    DO iPatch2 = 1,pGrid%nPatches
      pPatch2 => pRegion%patches(iPatch2)
    
      DO iPatch = iPatch2+1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        pPatch%nbMap(pPatch2%iPatchGlobal) = pPatch2%nbMap(pPatch%iPatchGlobal)
      END DO ! iPatch
    END DO ! iPatch2
  END IF ! pGrid%nPatches

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building patch-neighbor maps done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_BuildPatchNeighborMaps







! ******************************************************************************
!
! Purpose: Checking consistency of patch geometry with boundary conditions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_CheckPatchBcConsistency(pRegion)

  USE ModTools

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

  INTEGER :: errorFlag,iPatch,iPatchRelated
  REAL(RFREAL) :: angleSum,eqTol
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch,pPatchRelated

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckPatchBcConsistency',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking patch/bc consistency...' 
  END IF ! global%myProcid
 
  pGrid => pRegion%grid
 
  eqTol = 1.0E-12_RFREAL
 
! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    SELECT CASE ( pPatch%bcType ) 
      CASE ( BC_PERIODIC ) 
        iPatchRelated =  pPatch%iPatchRelated
        pPatchRelated => pRegion%patches(iPatchRelated)
        
        angleSum = pPatch%angleRelated+pPatchRelated%angleRelated
        
        IF ( FloatEqual(angleSum,0.0_RFREAL,eqTol) .EQV. .FALSE. ) THEN 
          CALL ErrorStop(global,ERR_PATCH_BC_INCONSISTENT,__LINE__)
        END IF ! FloatEqual
        
        IF ( pPatch%axisRelated /= pPatchRelated%axisRelated ) THEN 
          CALL ErrorStop(global,ERR_PATCH_BC_INCONSISTENT,__LINE__)
        END IF ! pPatch%axisRelated
        
        IF ( pPatch%axisRelated /= 1 .AND. & 
             pPatch%axisRelated /= 2 .AND. & 
             pPatch%axisRelated /= 3 ) THEN
          CALL ErrorStop(global,ERR_PATCH_BC_INCONSISTENT,__LINE__)   
        END IF ! pPatch%axisRelated
      CASE ( BC_SYMMETRY ) 
        IF ( pPatch%flatFlag .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_PATCH_BC_INCONSISTENT,__LINE__)
        END IF ! pPatch%flatFlag    
    END SELECT ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking patch/bc consistency done.'
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_CheckPatchBcConsistency





! *******************************************************************************
!
! Purpose: Compute global patch normal vectors and flatness flags.
!
! Description: None.
!
! Input:
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ComputePatchNormalsGlobal(regions)

  IMPLICIT NONE
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: xFlatFlag,yFlatFlag,zFlatFlag
  LOGICAL, DIMENSION(:), ALLOCATABLE :: flatFlag
  INTEGER :: errorFlag,ifl,iPatch,iReg
  REAL(RFREAL) :: eqTol,nx,nxMax,nxMin,ny,nyMax,nyMin,nz,nzMax,nzMin
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: globalVals,localVals
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pfnExt
  TYPE(t_global), POINTER :: global  
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion  

! ******************************************************************************
! Start, set pointers and variables
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_ComputePatchNormalsGlobal',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing global patch normal vectors...'                                                        
  END IF ! global%myProcid         

  eqTol = 1.0E-5_RFREAL

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,E15.8)') SOLVER_NAME,'Tolerance:',eqTol 
    WRITE(STDOUT,'(A,3X,A,2X,A,2X,A,2X,A)') SOLVER_NAME, 'Local','Global', &
                                            'Flat','Normal vector'
  END IF ! global%myProcid 

! ******************************************************************************
! Allocate and initialize temporary memory
! ******************************************************************************

  ALLOCATE(pfnExt(MIN_VAL:MAX_VAL,XCOORD:ZCOORD,global%nPatches), & 
           STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pfnExt')
  END IF ! global%error

  DO iPatch = 1,global%nPatches
    pfnExt(MIN_VAL,XCOORD,iPatch) =  HUGE(1.0_RFREAL)
    pfnExt(MIN_VAL,YCOORD,iPatch) =  HUGE(1.0_RFREAL)
    pfnExt(MIN_VAL,ZCOORD,iPatch) =  HUGE(1.0_RFREAL) 
    
    pfnExt(MAX_VAL,XCOORD,iPatch) = -HUGE(1.0_RFREAL)
    pfnExt(MAX_VAL,YCOORD,iPatch) = -HUGE(1.0_RFREAL)
    pfnExt(MAX_VAL,ZCOORD,iPatch) = -HUGE(1.0_RFREAL)  
  END DO ! iPatch     
     
  ALLOCATE(flatFlag(global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'flatFlag')
  END IF ! global%error     
     
! ******************************************************************************
! Compute local extrema of patch face-normal vectors
! ******************************************************************************
     
  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)
    pGrid   => pRegion%grid
       
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      nxMin = pfnExt(MIN_VAL,XCOORD,pPatch%iPatchGlobal)
      nyMin = pfnExt(MIN_VAL,YCOORD,pPatch%iPatchGlobal)
      nzMin = pfnExt(MIN_VAL,ZCOORD,pPatch%iPatchGlobal)

      nxMax = pfnExt(MAX_VAL,XCOORD,pPatch%iPatchGlobal)
      nyMax = pfnExt(MAX_VAL,YCOORD,pPatch%iPatchGlobal)
      nzMax = pfnExt(MAX_VAL,ZCOORD,pPatch%iPatchGlobal)

      DO ifl = 1,pPatch%nBFaces
        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)

        nxMax = MAX(nx,nxMax)
        nxMin = MIN(nx,nxMin)

        nyMax = MAX(ny,nyMax)
        nyMin = MIN(ny,nyMin)

        nzMax = MAX(nz,nzMax)
        nzMin = MIN(nz,nzMin)        
      END DO ! ifl
      
      pfnExt(MIN_VAL,XCOORD,pPatch%iPatchGlobal) = nxMin
      pfnExt(MIN_VAL,YCOORD,pPatch%iPatchGlobal) = nyMin
      pfnExt(MIN_VAL,ZCOORD,pPatch%iPatchGlobal) = nzMin

      pfnExt(MAX_VAL,XCOORD,pPatch%iPatchGlobal) = nxMax
      pfnExt(MAX_VAL,YCOORD,pPatch%iPatchGlobal) = nyMax
      pfnExt(MAX_VAL,ZCOORD,pPatch%iPatchGlobal) = nzMax
    END DO ! iPatch
  END DO ! iReg     
             
! ******************************************************************************
! Compute global extrema of patch face-normal vectors
! ******************************************************************************

  ALLOCATE(localVals(XCOORD:ZCOORD,global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals')
  END IF ! global%error 

  ALLOCATE(globalVals(XCOORD:ZCOORD,global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals')
  END IF ! global%error 

  DO iPatch = 1,global%nPatches
    localVals(XCOORD,iPatch) = pfnExt(MIN_VAL,XCOORD,iPatch)
    localVals(YCOORD,iPatch) = pfnExt(MIN_VAL,YCOORD,iPatch)
    localVals(ZCOORD,iPatch) = pfnExt(MIN_VAL,ZCOORD,iPatch)
  END DO ! iPatch
  
  CALL MPI_AllReduce(localVals,globalVals,(ZCOORD-XCOORD+1)*global%nPatches, &
                     MPI_RFREAL,MPI_MIN,global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  DO iPatch = 1,global%nPatches
    pfnExt(MIN_VAL,XCOORD,iPatch) = globalVals(XCOORD,iPatch)
    pfnExt(MIN_VAL,YCOORD,iPatch) = globalVals(YCOORD,iPatch)
    pfnExt(MIN_VAL,ZCOORD,iPatch) = globalVals(ZCOORD,iPatch)
        
    localVals(XCOORD,iPatch) = pfnExt(MAX_VAL,XCOORD,iPatch)
    localVals(YCOORD,iPatch) = pfnExt(MAX_VAL,YCOORD,iPatch)
    localVals(ZCOORD,iPatch) = pfnExt(MAX_VAL,ZCOORD,iPatch)
  END DO ! iPatch
   
  CALL MPI_AllReduce(localVals,globalVals,(ZCOORD-XCOORD+1)*global%nPatches, &
                     MPI_RFREAL,MPI_MAX,global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  DO iPatch = 1,global%nPatches
    pfnExt(MAX_VAL,XCOORD,iPatch) = globalVals(XCOORD,iPatch)
    pfnExt(MAX_VAL,YCOORD,iPatch) = globalVals(YCOORD,iPatch)
    pfnExt(MAX_VAL,ZCOORD,iPatch) = globalVals(ZCOORD,iPatch)
  END DO ! iPatch

  DEALLOCATE(localVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals')
  END IF ! global%error 

  DEALLOCATE(globalVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals')
  END IF ! global%error 

! ******************************************************************************
! Determine whether patches are flat
! ******************************************************************************
     
  DO iPatch = 1,global%nPatches
    nxMin = pfnExt(MIN_VAL,XCOORD,iPatch)
    nyMin = pfnExt(MIN_VAL,YCOORD,iPatch)
    nzMin = pfnExt(MIN_VAL,ZCOORD,iPatch)

    nxMax = pfnExt(MAX_VAL,XCOORD,iPatch)
    nyMax = pfnExt(MAX_VAL,YCOORD,iPatch)
    nzMax = pfnExt(MAX_VAL,ZCOORD,iPatch)

    CALL RFLU_SetPatchFlatFlags(global,nxMin,nxMax,nyMin,nyMax,nzMin,nzMax, &
                                eqTol,xFlatFlag,yFlatFlag,zFlatFlag)
    
    IF ( (xFlatFlag .EQV. .TRUE.) .AND. & 
         (yFlatFlag .EQV. .TRUE.) .AND. & 
         (zFlatFlag .EQV. .TRUE.) ) THEN 
      flatFlag(iPatch) = .TRUE.   
    ELSE 
      flatFlag(iPatch) = .FALSE.     
    END IF ! xFlatFlag   
  END DO ! iPatch

! ******************************************************************************
! Set patch normal and flatness flag
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)
    pGrid   => pRegion%grid
       
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
    
      IF ( flatFlag(pPatch%iPatchGlobal) .EQV. .TRUE. ) THEN 
        pPatch%flatFlag = .TRUE.

        pPatch%pn(XCOORD) = pfnExt(MAX_VAL,XCOORD,pPatch%iPatchGlobal)
        pPatch%pn(YCOORD) = pfnExt(MAX_VAL,YCOORD,pPatch%iPatchGlobal)
        pPatch%pn(ZCOORD) = pfnExt(MAX_VAL,ZCOORD,pPatch%iPatchGlobal)   
      ELSE 
        pPatch%flatFlag = .FALSE.

        pPatch%pn(XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%pn(YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%pn(ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)      
      END IF ! flatFlag
      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        IF ( pPatch%flatFlag .EQV. .TRUE. ) THEN  
          WRITE(STDOUT,'(A,2X,I4,4X,I4,6X,L1,3(2X,E13.6))') & 
            SOLVER_NAME,iPatch,pPatch%iPatchGlobal,pPatch%flatFlag, & 
            pPatch%pn(XCOORD:ZCOORD)
        ELSE 
          WRITE(STDOUT,'(A,2X,I4,4X,I4,6X,L1)') & 
            SOLVER_NAME,iPatch,pPatch%iPatchGlobal,pPatch%flatFlag       
        END IF ! pPatch%flatFlag
      END IF ! global%myProcid      
    END DO ! iPatch    
  END DO ! iReg

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(pfnExt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pfnExt')
  END IF ! global%error
       
  DEALLOCATE(flatFlag,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'flatFlag')
  END IF ! global%error     

! ******************************************************************************
! End  
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing global patch normal vectors done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePatchNormalsGlobal
  
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Compute global patch normal vectors and flatness flags.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ComputePatchNormalsLocal(pRegion)

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

  LOGICAL :: xFlatFlag,yFlatFlag,zFlatFlag
  INTEGER :: errorFlag,ifl,iPatch
  REAL(RFREAL) :: eqTol,nx,nxMax,nxMin,ny,nyMax,nyMin,nz,nzMax,nzMin
  TYPE(t_global), POINTER :: global  
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch  

! ******************************************************************************
! Start, set pointers and variables
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputePatchNormalsLocal',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing patch normal vectors...'                          
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal  
  END IF ! global%myProcid            
  
  eqTol = 1.0E-5_RFREAL  
  
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,E15.8)') SOLVER_NAME,'Tolerance:',eqTol          
    WRITE(STDOUT,'(A,3X,A,2X,A,2X,A,2X,A)') SOLVER_NAME,'Local','Global', & 
                                            'Flat','Normal vector'
  END IF ! global%myProcid 
  
  pGrid => pRegion%grid
  
! ******************************************************************************
! Compute patch normal vectors
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    nxMin =  HUGE(1.0_RFREAL)
    nyMin =  HUGE(1.0_RFREAL)
    nzMin =  HUGE(1.0_RFREAL)

    nxMax = -HUGE(1.0_RFREAL)
    nyMax = -HUGE(1.0_RFREAL)
    nzMax = -HUGE(1.0_RFREAL)

! ==============================================================================
!   Compute local extrema of patch face-normal vectors
! ==============================================================================

    DO ifl = 1,pPatch%nBFaces
      nx = pPatch%fn(XCOORD,ifl)
      ny = pPatch%fn(YCOORD,ifl)
      nz = pPatch%fn(ZCOORD,ifl)

      nxMax = MAX(nx,nxMax)
      nxMin = MIN(nx,nxMin)

      nyMax = MAX(ny,nyMax)
      nyMin = MIN(ny,nyMin)

      nzMax = MAX(nz,nzMax)
      nzMin = MIN(nz,nzMin)        
    END DO ! ifl

! ==============================================================================
!   Determine whether patches are flat
! ==============================================================================

    CALL RFLU_SetPatchFlatFlags(global,nxMin,nxMax,nyMin,nyMax,nzMin,nzMax, &
                                eqTol,xFlatFlag,yFlatFlag,zFlatFlag)

! ==============================================================================
!   Set patch normal and flatness flag
! ==============================================================================

    IF ( (xFlatFlag .EQV. .TRUE.) .AND. & 
         (yFlatFlag .EQV. .TRUE.) .AND. & 
         (zFlatFlag .EQV. .TRUE.) ) THEN 
      pPatch%flatFlag = .TRUE.
      
      pPatch%pn(XCOORD) = nxMax
      pPatch%pn(YCOORD) = nyMax
      pPatch%pn(ZCOORD) = nzMax      
    ELSE 
      pPatch%flatFlag = .FALSE.
      
      pPatch%pn(XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pPatch%pn(YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pPatch%pn(ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)      
    END IF ! xFlatFlag

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      IF ( pPatch%flatFlag .EQV. .TRUE. ) THEN  
        WRITE(STDOUT,'(A,2X,I4,4X,I4,6X,L1,3(2X,E13.6))') & 
          SOLVER_NAME,iPatch,pPatch%iPatchGlobal,pPatch%flatFlag, & 
          pPatch%pn(XCOORD:ZCOORD)
      ELSE 
        WRITE(STDOUT,'(A,2X,I4,4X,I4,6X,L1)') & 
          SOLVER_NAME,iPatch,pPatch%iPatchGlobal,pPatch%flatFlag       
      END IF ! pPatch%flatFlag
    END IF ! global%myProcid
  END DO ! iPatch
  
! ******************************************************************************
! End  
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing patch normal vectors done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePatchNormalsLocal
  
  
 



! ******************************************************************************
!
! Purpose: Create patch-neighbor maps.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_CreatePatchNeighborMaps(pRegion)

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

  INTEGER :: errorFlag,iPatch,iPatch2
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CreatePatchNeighborMaps',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating patch-neighbor maps...' 
  END IF ! global%myProcid
 
  pGrid => pRegion%grid
 
! ******************************************************************************
! Loop over patches and allocate memory
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
       
    ALLOCATE(pPatch%nbMap(global%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%nbMap')
    END IF ! global%error
    
    DO iPatch2 = 1,global%nPatches
      pPatch%nbMap(iPatch2) = .FALSE.
    END DO ! iPatch2
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating patch-neighbor maps done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_CreatePatchNeighborMaps
    
  







! ******************************************************************************
!
! Purpose: Destroy patch-neighbor maps.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_DestroyPatchNeighborMaps(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyPatchNeighborMaps',&
  'RFLU_ModPatchUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying patch-neighbor maps...' 
  END IF ! global%myProcid
 
  pGrid => pRegion%grid
 
! ******************************************************************************
! Loop over patches and deallocate memory
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
       
    DEALLOCATE(pPatch%nbMap,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%nbMap')
    END IF ! global%error
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying patch-neighbor maps done.' 
  END IF ! global%myProcid
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_DestroyPatchNeighborMaps 
 



! ******************************************************************************
!
! Purpose: Determine whether patches are flat by checking normal extrema.
!
! Description: Checking whether extrema have same sign (if yes, and extrema are 
!   equal, have flat patch; if not, and extrema are close to zero, also have 
!   flat patch)
!
! Input: 
!   global      Pointer to global data
!   nxMin       Minimum of x-component of normal vector
!   nxMax       Maximum of x-component of normal vector
!   nyMin       Minimum of y-component of normal vector
!   nyMax       Maximum of y-component of normal vector
!   nzMin       Minimum of z-component of normal vector
!   nzMax       Maximum of z-component of normal vector
!   eqTol       Equality tolerance
!
! Output: 
!   xFlatFlag   Flatness flag for x-component
!   yFlatFlag   Flatness flag for y-component
!   zFlatFlag   Flatness flag for z-component
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_SetPatchFlatFlags(global,nxMin,nxMax,nyMin,nyMax,nzMin,nzMax, &
                                  eqTol,xFlatFlag,yFlatFlag,zFlatFlag)

  USE ModTools

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(OUT) :: xFlatFlag,yFlatFlag,zFlatFlag
  REAL(RFREAL), INTENT(IN) :: eqTol,nxMax,nxMin,nyMax,nyMin,nzMax,nzMin
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_SetPatchFlatFlags',&
  'RFLU_ModPatchUtils.F90')

! ******************************************************************************
! Determine flatness flags
! ******************************************************************************

  IF ( (SIGN(1.0_RFREAL,nxMin) == SIGN(1.0_RFREAL,nxMax)) ) THEN 
    IF ( FloatEqual(nxMin,nxMax,eqTol) .EQV. .TRUE. ) THEN 
      xFlatFlag = .TRUE.
    ELSE 
      xFlatFlag = .FALSE.
    END IF ! FloatEqual
  ELSE
    IF ( (FloatEqual(ABS(nxMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
         (FloatEqual(ABS(nxMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
      xFlatFlag = .TRUE.      
    ELSE
      xFlatFlag = .FALSE.
    END IF ! FloatEqual
  END IF ! FloatEqual

  IF ( (SIGN(1.0_RFREAL,nyMin) == SIGN(1.0_RFREAL,nyMax)) ) THEN 
    IF ( FloatEqual(nyMin,nyMax,eqTol) .EQV. .TRUE. ) THEN 
      yFlatFlag = .TRUE.
    ELSE 
      yFlatFlag = .FALSE.
    END IF ! FloatEqual
  ELSE
    IF ( (FloatEqual(ABS(nyMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
         (FloatEqual(ABS(nyMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
      yFlatFlag = .TRUE.      
    ELSE
      yFlatFlag = .FALSE.
    END IF ! FloatEqual
  END IF ! FloatEqual

  IF ( (SIGN(1.0_RFREAL,nzMin) == SIGN(1.0_RFREAL,nzMax)) ) THEN     
    IF ( FloatEqual(nzMin,nzMax,eqTol) .EQV. .TRUE. ) THEN 
      zFlatFlag = .TRUE.
    ELSE 
      zFlatFlag = .FALSE.
    END IF ! FloatEqual
  ELSE    
    IF ( (FloatEqual(ABS(nzMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
         (FloatEqual(ABS(nzMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
      zFlatFlag = .TRUE.      
    ELSE
      zFlatFlag = .FALSE.
    END IF ! FloatEqual
  END IF ! FloatEqual

! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_SetPatchFlatFlags






! ******************************************************************************
!
! Purpose: Get patch normal direction.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   pPatch      Pointer to patch data
!
! Output: 
!   pnDir       Patch normal direction (equal to coordinate index if normal
!               aligned with coordinate direction, otherwise set to 
!               CRAZY_VALUE_INT)
!   pnDirFlag   Patch normal direction flag (TRUE if normal aligned with 
!               coordinate axis, FALSE otherwise)
!
! Notes: None.
!
! ******************************************************************************
  
SUBROUTINE RFLU_GetPatchNormalDirection(global,pPatch,pnDir,pnDirFlag)

  USE ModTools

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(OUT) :: pnDirFlag
  INTEGER, INTENT(OUT) :: pnDir
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: pnDirTemp(1)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_GetPatchNormalDirection',&
  'RFLU_ModPatchUtils.F90')

! ******************************************************************************
! 
! ******************************************************************************

  IF ( pPatch%flatFlag .EQV. .FALSE. ) THEN 
    pnDir     = CRAZY_VALUE_INT
    pnDirFlag = .FALSE.
  ELSE 
    pnDirTemp = MAXLOC(ABS(pPatch%pn(XCOORD:ZCOORD)))       
    pnDir     = NINT(pnDirTemp(1))
    pnDirFlag = FloatEqual(ABS(pPatch%pn(pnDir)),1.0_RFREAL,1.0E-6_RFREAL)
  END IF ! pPatch%flatFlag

! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global) 
 
END SUBROUTINE RFLU_GetPatchNormalDirection





END MODULE RFLU_ModPatchUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPatchUtils.F90,v $
! Revision 1.7  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/12/21 12:21:43  haselbac
! Bug fix: Wrong operation argument to MPI call
!
! Revision 1.4  2006/06/14 20:11:24  mparmar
! Bug fix: need to use iPatchGlobal instead of iPatch
!
! Revision 1.3  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.2  2006/04/07 14:50:16  haselbac
! Added RFLU_GetPatchNormalDirection
!
! Revision 1.1  2006/03/25 21:38:54  haselbac
! Initial revision
!
! ******************************************************************************














