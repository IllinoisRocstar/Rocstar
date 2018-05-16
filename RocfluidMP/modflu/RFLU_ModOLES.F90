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
!*******************************************************************************
!
! Purpose: Suite of routines for optimal LES computations.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: RFLU_ModOLES.F90,v 1.10 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModOLES

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModTools, ONLY: MakeNonZero
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModSortSearch
  USE ModTools
  USE ModMPI
  
  USE RFLU_ModOctree

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_CreateStencilsWeightsOLES, & 
            RFLU_CreateIntegralsOLES, & 
            RFLU_BuildStencilsOLES, & 
            RFLU_FindPrototypeFacesOLES, & 
            RFLU_ComputeGeometricTermsOLES, & 
            RFLU_BuildSymmetryMapsOLES, & 
            RFLU_EnforceSymmetryOLES, & 
            RFLU_DefineCorrelation220, & 
            RFLU_DefineCorrelation221, & 
            RFLU_DefineCorrelation32, & 
            RFLU_DefineCorrelation430, & 
            RFLU_DefineCorrelation431, & 
            RFLU_DefineCorrelation432, & 
            RFLU_DefineCorrelation540, & 
            RFLU_DefineCorrelation541, & 
            RFLU_DefineCorrelation542, & 
            RFLU_MapK2IJ, & 
            RFLU_MapL2IJK, & 
            RFLU_MapM2IJKL, & 
            RFLU_GetI1PosOLES, & 
            RFLU_GetI4PosOLES, & 
            RFLU_GetLPosOLES, & 
            RFLU_GetLPosInvOLES, & 
            RFLU_GetQPosOLES, & 
            RFLU_GetQPosInvOLES, & 
            RFLU_DestroyStencilsWeightsOLES, & 
            RFLU_AllocateDCUHREArrays, & 
            RFLU_DeallocateDCUHREArrays, &
            RFLU_SetMapFunNZ2FunCorr22, & 
            RFLU_SetMapFunNZ2FunCorr32, & 
            RFLU_SetMapFunNZ2FunCorr43, & 
            RFLU_SetMapFunNZ2FunCorr44
            
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModOLES.F90,v $ $Revision: 1.10 $'        
  TYPE(t_grid), POINTER :: pGrid     
    
! ==============================================================================
! Variables and parameters used by DCUHRE subroutine
! ==============================================================================      
       
  INTEGER, PUBLIC :: maxCalls,nDim,nEval,workArraySize,workArraySizeNew    
  INTEGER, PARAMETER, PUBLIC :: DCUHRE_LOOP_LIMIT = 5,        &  
                                MAX_CALLS_FACTOR  = 5,        &
                                MAX_CALLS_LIMIT   = 50000000, & 
                                MAX_CALLS_START   = 100000,   &
                                MIN_CALLS         = 0
  INTEGER, PUBLIC :: dummy(1)
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: symMapI45OLES  
    
  REAL(RFREAL), PUBLIC :: errAbsReq,errRelReq
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE, PUBLIC :: lowLim,errAbsEst, &
                                                     uppLim,integral, &
                                                     integralNZ,workArray                  
              
! ==============================================================================
! Variables and parameters used in integration routines
! ==============================================================================       
       
  INTEGER, PUBLIC :: nzLoc,nzSgn 
  INTEGER, DIMENSION(3,3), PARAMETER, PUBLIC :: mapSurf2Vol2 = & 
    RESHAPE((/6,4,4,4,6,5,5,5,6/),(/3,3/))
  INTEGER, DIMENSION(3,6), PARAMETER, PUBLIC :: mapSurf2Vol3 = & 
    RESHAPE((/9,7,7,7,9,8,8,8,9,9,7,7,7,9,8,8,8,9/),(/3,6/))
  INTEGER, DIMENSION(3,3), PARAMETER, PUBLIC :: kd = & ! Kronecker delta
    RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))     
        
  INTEGER, PUBLIC :: mapFunNZ2FunCorr22(9),mapFunNZ2FunCorr32(9), & 
                     mapFunNZ2FunCorr43(27),mapFunNZ2FunCorr44(81)  
    
  REAL(RFREAL), PUBLIC :: nzMag,nzVal 
  REAL(RFREAL), PARAMETER, PUBLIC :: CONST_KOLMOGOROV = 2.0_RFREAL  
       

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
! ******************************************************************************
!   Create geometry
! ******************************************************************************

    SUBROUTINE RFLU_CreateStencilsWeightsOLES(pRegion)
  
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables

      INTEGER :: errorFlag,nCells,nFaces
      TYPE(t_grid), POINTER :: pGrid           
      TYPE(t_global), POINTER :: global      
            
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_CreateStencilsWeightsOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating optimal LES '// & 
                                 'stencils, weights, and mappings...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer
! ==============================================================================

      pGrid => pRegion%grid

! ==============================================================================
!     Allocate memory for face stencils, cell limits
! ==============================================================================

! TEMPORARY: hard-code stencil
      nCells = 2
! END TEMPORARY

      ALLOCATE(pGrid%fsOLES(nCells,pGrid%nFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%fsOLES')
      END IF ! global%error

      pGrid%fsOLES(:,:) = 0

      ALLOCATE(pGrid%intLimOLES(INT_LIM_LOW:INT_LIM_UPP,XCOORD:ZCOORD, & 
               pGrid%nCellsTot),STAT=errorFlag)
      global%error = errorFlag            
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%intLimOLES')
      END IF ! global%error

      pGrid%intLimOLES(INT_LIM_LOW,:,:) =  HUGE(1.0_RFREAL)
      pGrid%intLimOLES(INT_LIM_UPP,:,:) = -HUGE(1.0_RFREAL)

      ALLOCATE(pGrid%rhoOLES(pGrid%nFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%intLimOLES')
      END IF ! global%error

      pGrid%rhoOLES(:) = HUGE(1.0_RFREAL)

      ALLOCATE(pGrid%fp2fOLES(3),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%fp2fOLES')
      END IF ! global%error

      ALLOCATE(pGrid%f2fpOLES(pGrid%nFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%f2fpOLES')
      END IF ! global%error

! ==============================================================================
!     Weights, NOTE store only for prototype faces
! ==============================================================================

      nFaces = 3

      ALLOCATE(pGrid%wtLinOLES(3,3,nCells,nFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%wtLinOLES')
      END IF ! global%error      

      ALLOCATE(pGrid%wtQuadOLES(3,3,3,nCells,nCells,nFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%wQuadOLES')
      END IF ! global%error   

      pGrid%wtLinOLES(:,:,:,:)      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%wtQuadOLES(:,:,:,:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ==============================================================================
!     Symmetry maps for integrals 
! ==============================================================================

      ALLOCATE(symMapI45OLES(9*nCells*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'symMapI45OLES')      
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating optimal LES '// & 
                                 'stencils, weights, and mappings done.'
      END IF ! global%verbLevel
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_CreateStencilsWeightsOLES



! ******************************************************************************
!   Create integrals. NOTE initialize to crazy values to be able to check 
!   whether values are entered correctly.
! ******************************************************************************

    SUBROUTINE RFLU_CreateIntegralsOLES(pRegion)
  
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables

      INTEGER :: errorFlag,nCells,nCols,nFaces,nRows
      TYPE(t_grid), POINTER :: pGrid           
      TYPE(t_global), POINTER :: global             
            
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_CreateIntegralsOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating optimal LES '// & 
                                 'integrals...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer
! ==============================================================================

      pGrid => pRegion%grid

! ==============================================================================
!     Allocate memory for integrals, NOTE store only for prototype faces
! ==============================================================================

      nFaces = 3
      nCells = SIZE(pGrid%fsOLES,1)

! ------------------------------------------------------------------------------
!     Integral 1 
! ------------------------------------------------------------------------------

      ALLOCATE(pGrid%int1OLES(3,nFaces,3*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int1OLES')
      END IF ! global%error

      pGrid%int1OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ------------------------------------------------------------------------------
!     Integral 2 
! ------------------------------------------------------------------------------

      ALLOCATE(pGrid%int20OLES(nFaces,3*nCells,3*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int20OLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%int21OLES(nFaces,3*nCells,3*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int21OLES')
      END IF ! global%error     

      pGrid%int20OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int21OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)


! ------------------------------------------------------------------------------
!     Integral 3 
! ------------------------------------------------------------------------------

      ALLOCATE(pGrid%int31OLES(nFaces,3*nCells,9*nCells*nCells), & 
               STAT=errorFlag)
      global%error = errorFlag            
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int31OLES')
      END IF ! global%error

      ALLOCATE(pGrid%int32OLES(nFaces,9*nCells*nCells,3*nCells), & 
               STAT=errorFlag)           
      global%error = errorFlag           
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int32OLES')
      END IF ! global%error
   
      pGrid%int31OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int32OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ------------------------------------------------------------------------------
!     Integral 4 
! ------------------------------------------------------------------------------

      ALLOCATE(pGrid%int40OLES(3,nFaces,9*nCells*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int40OLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%int41OLES(3,nFaces,9*nCells*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int41OLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%int42OLES(3,nFaces,9*nCells*nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int42OLES')
      END IF ! global%error           

      pGrid%int40OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int41OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int42OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)            

! ------------------------------------------------------------------------------
!     Integral 5
! ------------------------------------------------------------------------------

      ALLOCATE(pGrid%int50OLES(nFaces,9*nCells*nCells,9*nCells*nCells), & 
               STAT=errorFlag)
      global%error = errorFlag           
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int50OLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%int51OLES(nFaces,9*nCells*nCells,9*nCells*nCells), & 
               STAT=errorFlag)
      global%error = errorFlag            
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int51OLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%int52OLES(nFaces,9*nCells*nCells,9*nCells*nCells), & 
               STAT=errorFlag)
      global%error = errorFlag            
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%int52OLES')
      END IF ! global%error

      pGrid%int50OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int51OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pGrid%int52OLES(:,:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)            

! ------------------------------------------------------------------------------
!     LHS and RHS of linear system
! ------------------------------------------------------------------------------

      nRows = 3*nCells*(1 + 3*nCells)
      nCols = nRows

      ALLOCATE(pGrid%lhsOLES(nRows,nCols),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%lhsOLES')
      END IF ! global%error

      ALLOCATE(pGrid%lhsInvOLES(nCols,nRows),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%lhsInvOLES')
      END IF ! global%error
      
      ALLOCATE(pGrid%rhsOLES(nRows),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%rhsOLES')
      END IF ! global%error
      
      pGrid%lhsOLES(:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)      
      pGrid%rhsOLES(:)   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ==============================================================================
!     End
! ==============================================================================

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating optimal LES '// & 
                                 'integrals done.'
      END IF ! global%verbLevel
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_CreateIntegralsOLES




! ******************************************************************************
!   Build stencil
! ******************************************************************************

    SUBROUTINE RFLU_BuildStencilsOLES(pRegion)
   
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables

      INTEGER :: errorFlag,ic,ifc,is,loc,locOffset1,locOffset2,stencilSize, & 
                 stencilSizeTemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: stencilTemp
      REAL(RFREAL) :: delFrac,dotp,vecm,xDel,xMax,xMin,yDel,yMax,yMin, & 
                      zDel,zMax,zMin
      REAL(RFREAL) :: dummy(1),vec(3)
      REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: dist  
      TYPE(t_grid), POINTER :: pGrid              
      TYPE(t_global), POINTER :: global     
    
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_BuildStencilsOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building optimal LES stencils...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer and constants
! ==============================================================================

      pGrid => pRegion%grid
      
      delFrac = 0.01_RFREAL      
      
! ==============================================================================
!     Build face stencils
! ==============================================================================      
      
      stencilSize = SIZE(pGrid%fsOLES,1)           
                      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Building face stencil support...'
        WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Stencil width:',stencilSize
      END IF ! global%verbLevel               
           
! ------------------------------------------------------------------------------
!     Stencil only contains cells which straddle face, Octree not necessary
! ------------------------------------------------------------------------------

      IF ( stencilSize == 2 ) THEN ! only need cells which straddle face
        DO ifc = 1,pGrid%nFaces
          pGrid%fsOLES(1,ifc) = pGrid%f2c(1,ifc)
          pGrid%fsOLES(2,ifc) = pGrid%f2c(2,ifc)
        END DO ! ifc

! ------------------------------------------------------------------------------
!     Stencil contains more cells, use Octree
! ------------------------------------------------------------------------------

      ELSE IF ( stencilSize > 2 ) THEN 
        IF ( stencilSize == 4 ) THEN 
          stencilSizeTemp = 20 ! to make sure to capture distance-two cells
        ELSE IF ( stencilSize == 6 ) THEN 
          stencilSizeTemp = 102 ! to make sure to capture distance-three cells
        ELSE IF ( stencilSize == 8 ) THEN 
          stencilSizeTemp = 296 ! to make sure to capture distance-four cells
        ELSE   
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! stencilSize

! ----- Allocate memory -------------------------------------------------------

        ALLOCATE(stencilTemp(stencilSizeTemp),STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'stencilTemp')
        END IF ! global%error

        ALLOCATE(dist(stencilSizeTemp),STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dist')
        END IF ! global%error

! ----- Determine bounding box ------------------------------------------------

        xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
        xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
        yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
        yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
        zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))
        zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))

        xDel = xMax - xMin 
        yDel = yMax - yMin
        zDel = zMax - zMin

        xMin = xMin - delFrac*xDel
        xMax = xMax + delFrac*xDel 
        yMin = yMin - delFrac*yDel
        yMax = yMax + delFrac*yDel 
        zMin = zMin - delFrac*zDel
        zMax = zMax + delFrac*zDel 

! ----- Create and build Octree -----------------------------------------------

        CALL RFLU_CreateOctree(global,pGrid%nCellsTot)
        CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCellsTot), & 
                              pGrid%cofg(YCOORD,1:pGrid%nCellsTot), &
                              pGrid%cofg(ZCOORD,1:pGrid%nCellsTot), & 
                              xMin,xMax,yMin,yMax,zMin,zMax)

! ----- Loop over faces -------------------------------------------------------

        DO ifc = 1,pGrid%nFaces
          CALL RFLU_QueryOctree(pGrid%fc(XCOORD,ifc),pGrid%fc(YCOORD,ifc), & 
                                pGrid%fc(ZCOORD,ifc),stencilSizeTemp, & 
                                stencilTemp)

! ------- Compute distance (signed) normal to face. At present, this only 
!         accepts cells if they are located in the normal direction of the face,
!         so you can only generate 1x1x<nCells> stencils with this method. Not 
!         hard to get other stencils, but hard to control them on non-uniform 
!         grids. Rejected faces get a distance value of HUGE(1.0_RFREAL), which
!         is used further below in checking that enough cells were found.

          dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))
          loc   = dummy(1)    

          DO is = 1,stencilSizeTemp
            ic = stencilTemp(is)      

            vec(:) = pGrid%cofg(:,ic) - pGrid%fc(:,ifc)
            vecm   = SQRT(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
            vec(:) = vec(:)/vecm

            dotp = vec(1)*pGrid%fn(1,ifc) & 
                 + vec(2)*pGrid%fn(2,ifc) & 
                 + vec(3)*pGrid%fn(3,ifc) 

            IF ( FloatEqual(ABS(dotp),1.0_RFREAL) .EQV. .TRUE. ) THEN        
              dist(is) = pGrid%cofg(loc,ic) - pGrid%fc(loc,ifc)
            ELSE 
              dist(is) = HUGE(1.0_RFREAL)
            END IF ! dist  
          END DO ! is

! ------- Sort cells according to increasing distance

          CALL QuickSortRFREALInteger(dist(1:stencilSizeTemp), & 
                                      stencilTemp(1:stencilSizeTemp), & 
                                      stencilSizeTemp)          
                                    
! ------- Make sure that stencil contains cells which straddle face, and that
!         the two cells are in the correct location (very strict check). 

          IF ( pGrid%fn(loc,ifc) > 0.0_RFREAL ) THEN 
            locOffset1 = 0
            locOffset2 = 1
          ELSE 
            locOffset1 = 1
            locOffset2 = 0
          END IF ! pGrid%fn

          DO ic = 1,2            
            IF ( ic == 1 ) THEN ! NOTE integer division
              loc = stencilSize/2 + locOffset1
            ELSE 
              loc = stencilSize/2 + locOffset2
            END IF ! ic
            
            IF ( stencilTemp(loc) /= pGrid%f2c(ic,ifc) ) THEN 
              CALL ErrorStop(global,ERR_OLES_STENCIL,__LINE__)              
            END IF ! stencilTemp                         
          END DO ! ic

! ------- Make sure that the first <stencilSize> cells are in line with the face,
!         in other words, do not have dist(is) = HUGE(1.0_RFREAL)

          DO is = 1,stencilSize
            IF ( FloatEqual(dist(is),HUGE(1.0_RFREAL)) .EQV. .TRUE. ) THEN 
              CALL ErrorStop(global,ERR_OLES_STENCIL,__LINE__)
            END IF ! dist
          END DO ! is

! ------- Copy first <stencilSize> cells into face stencil array

          DO is = 1,stencilSize          
            pGrid%fsOLES(is,ifc) = stencilTemp(is)
          END DO ! is
        END DO ! ifc            
      
! ----- Destroy Octree and deallocate memory ----------------------------------------                             
                              
        CALL RFLU_DestroyOctree(global) 

        DEALLOCATE(stencilTemp,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'stencilTemp')
        END IF ! global%error 

        DEALLOCATE(dist,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dist')
        END IF ! global%error    
      END IF ! stencilSize
      

#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES face stencils'
      IF ( stencilSize == 2 ) THEN 
        DO ifc = 1,pGrid%nFaces
          WRITE(STDOUT,'(A,3(1X,I6))') SOLVER_NAME,ifc,pGrid%fsOLES(1:2,ifc)
        END DO ! ifc
      ELSE IF ( stencilSize == 4 ) THEN 
        DO ifc = 1,pGrid%nFaces
          WRITE(STDOUT,'(A,5(1X,I6))') SOLVER_NAME,ifc,pGrid%fsOLES(1:4,ifc)
        END DO ! ifc        
      END IF ! stencilSize
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME         
#endif   
       
! ==============================================================================
!     End
! ==============================================================================
  
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building optimal LES '// &
                                 'stencils done.'
      END IF ! global%verbLevel   
  
      CALL DeregisterFunction(global)     
   
    END SUBROUTINE RFLU_BuildStencilsOLES




! ******************************************************************************
!   Find prototypical faces: to reduce computational cost 
! ******************************************************************************

    SUBROUTINE RFLU_FindPrototypeFacesOLES(pRegion)
   
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables

      INTEGER, PARAMETER :: LOCSAVE_INIT = -1 ! Must be zero or negative
      INTEGER :: ifc,ifcp,loc,locSaveCntr
      INTEGER :: dummy(1),locSave(3)
      TYPE(t_grid), POINTER :: pGrid              
      TYPE(t_global), POINTER :: global     
    
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_FindPrototypeFacesOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finding optimal LES '// & 
                                 'prototype faces...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer
! ==============================================================================

      pGrid => pRegion%grid
   
! ==============================================================================
!     Initialize
! ==============================================================================   
   
      locSaveCntr  = 0 
      locSave(1:3) = LOCSAVE_INIT
   
! ==============================================================================
!     Find prototype faces
! ==============================================================================      
      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Finding prototype faces...'
      END IF ! global%verbLevel      
      
      DO ifc = 1,pGrid%nFaces
        dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))
        loc   = dummy(1)        
      
        IF ( locSave(loc) == LOCSAVE_INIT ) THEN 
          locSave(loc) = 1
          locSaveCntr  = locSaveCntr + 1
          pGrid%fp2fOLES(loc) = ifc
        END IF ! locSave
      
        IF ( locSaveCntr == 3 ) THEN 
          EXIT
        END IF ! locSaveCntr
      END DO ! ifc
      
#ifdef CHECK_DATASTRUCT  
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Prototype faces:'
      DO ifcp = 1,3
        ifc = pGrid%fp2fOLES(ifcp)
        WRITE(STDOUT,'(A,2(1X,I6),3(1X,E18.9))') SOLVER_NAME,ifcp,ifc, & 
                                                 pGrid%fn(1:3,ifc)
      END DO ! ifcp  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME                 
#endif

! ==============================================================================
!     Map other faces onto prototype faces
! ==============================================================================

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mapping faces to '// & 
                                 'prototype faces...'
      END IF ! global%verbLevel    

      DO ifc = 1,pGrid%nFaces
        dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))
        loc   = dummy(1)         
      
        pGrid%f2fpOLES(ifc) = loc
      END DO ! ifc

#ifdef CHECK_DATASTRUCT  
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Prototype faces:'
      DO ifc = 1,pGrid%nFaces
        ifcp = pGrid%f2fpOLES(ifc)
        WRITE(STDOUT,'(A,2(1X,I6),3(1X,E18.9))') SOLVER_NAME,ifc,ifcp, & 
                                                 pGrid%fn(1:3,ifc)
      END DO ! ifc
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME                 
#endif

! ==============================================================================
!     End
! ==============================================================================
  
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finding optimal LES '// & 
                                 'prototype faces done.'
      END IF ! global%verbLevel   
  
      CALL DeregisterFunction(global)     
   
    END SUBROUTINE RFLU_FindPrototypeFacesOLES








! ******************************************************************************
!   Build stencil
! ******************************************************************************

   SUBROUTINE RFLU_ComputeGeometricTermsOLES(pRegion)
   
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables

      INTEGER :: ic,icl,icg,ifc,ip,c1,c2,loc,nCells 
      INTEGER :: dummy(1)
      REAL(RFREAL) :: term
      REAL(RFREAL) :: fc(XCOORD:ZCOORD)
      REAL(RFREAL) :: xyz(INT_LIM_LOW:INT_LIM_UPP,XCOORD:ZCOORD)
      TYPE(t_grid), POINTER :: pGrid              
      TYPE(t_patch), POINTER :: pPatch              
      TYPE(t_global), POINTER :: global             
            
#ifdef CHECK_DATASTRUCT
      INTEGER :: i,j
#endif      
            
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_ComputeGeometricTermsOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                                 'geometric terms...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer
! ==============================================================================

      pGrid => pRegion%grid

! ==============================================================================
!     Compute integration limits for each cell
! ============================================================================== 
   
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing cell integration '// & 
                                 'limits...'
        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Interior faces...'
        END IF ! global%verbLevel
      END IF ! global%verbLevel   
   
      DO ifc = 1,pGrid%nFacesTot ! NOTE loop over ALL faces
        c1 = pGrid%f2c(1,ifc)
        c2 = pGrid%f2c(2,ifc)
                
        dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))
        loc   = dummy(1)        
        
        fc(:) = ABS(pGrid%fn(loc,ifc))*pGrid%fc(:,ifc)
                
        IF ( c2 /= F2C_INIT ) THEN 
          IF ( pGrid%fn(loc,ifc) > 0.0_RFREAL ) THEN 
            pGrid%intLimOLES(INT_LIM_UPP,loc,c1) = fc(loc)
            pGrid%intLimOLES(INT_LIM_LOW,loc,c2) = fc(loc)
          ELSE   
            pGrid%intLimOLES(INT_LIM_LOW,loc,c1) = fc(loc)
            pGrid%intLimOLES(INT_LIM_UPP,loc,c2) = fc(loc)          
          END IF ! pGrid%intLimOLES
        ELSE  
          IF ( pGrid%fn(loc,ifc) > 0.0_RFREAL ) THEN 
            pGrid%intLimOLES(INT_LIM_UPP,loc,c1) = fc(loc)
          ELSE   
            pGrid%intLimOLES(INT_LIM_LOW,loc,c1) = fc(loc)          
          END IF ! pGrid%intLimOLES       
        END IF ! c2
      END DO ! ifc         
      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Patches...'
      END IF ! global%verbLevel      
      
      DO ip = 1,pGrid%nPatches
        pPatch => pRegion%patches(ip)

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,7X,A,I4)') SOLVER_NAME,'Patch: ',ip
        END IF ! global%verbLevel  

        DO ifc = 1,pPatch%nBFaces
          c1 = pPatch%bf2c(ifc)
        
          dummy = MAXLOC(ABS(pPatch%fn(1:3,ifc)))
          loc   = dummy(1)      
        
          fc(:) = ABS(pPatch%fn(loc,ifc))*pPatch%fc(:,ifc)
        
          IF ( pPatch%fn(loc,ifc) > 0.0_RFREAL ) THEN 
            pGrid%intLimOLES(INT_LIM_UPP,loc,c1) = fc(loc)
          ELSE   
            pGrid%intLimOLES(INT_LIM_LOW,loc,c1) = fc(loc)          
          END IF ! pGrid%intLimOLES       
        END DO ! ifc      
      END DO ! ip
      
      
#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES cell integration limits'
      DO ic = 1,pGrid%nCellsTot
        WRITE(STDOUT,'(A,1X,I6,6(1X,E18.9))') SOLVER_NAME,ic, & 
          ((pGrid%intLimOLES(i,j,ic),i=INT_LIM_LOW,INT_LIM_UPP),j=XCOORD,ZCOORD)
      END DO ! ic
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME         
#endif        
   
! ==============================================================================
!     Compute radius of smallest sphere contained in cells of stencil
! ============================================================================== 
   
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing smallest sphere '// & 
                                 'for stencils...'
      END IF ! global%verbLevel      
   
      DO ifc = 1,pGrid%nFaces    
        nCells = SIZE(pGrid%fsOLES,1)
      
        xyz(INT_LIM_LOW,XCOORD:ZCOORD) =  HUGE(1.0_RFREAL)      
        xyz(INT_LIM_UPP,XCOORD:ZCOORD) = -HUGE(1.0_RFREAL)
         
        DO icl = 1,nCells              
          icg = pGrid%fsOLES(icl,ifc)
          
          xyz(INT_LIM_LOW,1:3) = MIN(pGrid%intLimOLES(INT_LIM_LOW,1:3,icg), & 
                                     xyz(INT_LIM_LOW,1:3))
          xyz(INT_LIM_UPP,1:3) = MAX(pGrid%intLimOLES(INT_LIM_UPP,1:3,icg), & 
                                     xyz(INT_LIM_UPP,1:3))
        END DO ! ic        
 
! BEGIN ALTERNATIVE 1: Wrong, represents sphere containED in cells
        DO ic = 1,3
          pGrid%rhoOLES(ifc) = MIN(pGrid%rhoOLES(ifc), & 
                                   xyz(INT_LIM_UPP,ic)-xyz(INT_LIM_LOW,ic))
        END DO ! ic
! END ALTERNATIVE 1  
  
! BEGIN ALERNATIVE 2: Correct, represents sphere containing cells 
!                     Changes integrals, but not weights  
!        term = -HUGE(1.0_RFREAL)
!  
!        DO ic = 1,3
!          term = MAX(term,xyz(INT_LIM_UPP,ic)-xyz(INT_LIM_LOW,ic))
!        END DO ! ic
!        
!        pGrid%rhoOLES(ifc) = MIN(pGrid%rhoOLES(ifc),term)
! END ALTERNATIVE 2
      END DO ! ifc
              
#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES smallest stencil sphere'
      DO ifc = 1,pGrid%nFaces      
        WRITE(STDOUT,'(A,1X,I6,1X,E18.9)') SOLVER_NAME,ifc,pGrid%rhoOLES(ifc)
      END DO ! ifc
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME         
#endif               
      
! ==============================================================================
!     Compute filter width - assume constant grid spacing for now
! ==============================================================================      
      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing filter width...'
      END IF ! global%verbLevel      
            
      pGrid%deltaOLES = pGrid%vol(1)**(1.0_RFREAL/3.0_RFREAL)
         
#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES filter width'
      WRITE(STDOUT,'(A,1X,E18.9)') SOLVER_NAME,pGrid%deltaOLES
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME         
#endif         
      
! ==============================================================================
!     End
! ==============================================================================
  
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                                 'geometric terms done.'
      END IF ! global%verbLevel   
  
      CALL DeregisterFunction(global)     
   
    END SUBROUTINE RFLU_ComputeGeometricTermsOLES




! ******************************************************************************
!   Build symmetry maps for integrals
! ******************************************************************************

    SUBROUTINE RFLU_BuildSymmetryMapsOLES(pRegion)
    
      IMPLICIT NONE
      
! --- parameters      
      
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables    
    
      INTEGER :: b,g,j,k,nCells,pos,posSym
      TYPE(t_global), POINTER :: global     
    
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_BuildSymmetryMapsOLES',&
  'RFLU_ModOLES.F90')    
    
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A)') SOLVER_NAME            
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building optimal LES '// & 
                                 'symmetry maps...'
      END IF ! global%verbLevel    
    
! ==============================================================================
!     Build symmetry maps
! ==============================================================================     
 
      nCells = SIZE(pRegion%grid%fsOLES,1)
    
      DO b = 1,nCells
        DO g = 1,nCells  
          DO j = 1,3
            DO k = 1,3
              pos    = RFLU_GetQPosOLES(j,k,b,g,nCells)
              posSym = RFLU_GetQPosOLES(k,j,g,b,nCells) 

              IF ( posSym < pos ) THEN 
                symMapI45OLES(pos) = posSym
              ELSE 
                symMapI45OLES(pos) = 0        
              END IF ! pos2
            END DO ! k
          END DO ! j
        END DO ! g
      END DO ! b
    
! ==============================================================================
!     End
! ==============================================================================
  
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building optimal LES '// & 
                                 'symmetry maps done.'
      END IF ! global%verbLevel   
  
      CALL DeregisterFunction(global)      
    
    END SUBROUTINE RFLU_BuildSymmetryMapsOLES
    
    
    
    
! ******************************************************************************
!   Enforce symmetry: NOTE will not be completely symmetric!
! ******************************************************************************

    SUBROUTINE RFLU_EnforceSymmetryOLES(pRegion)
   
      IMPLICIT NONE
      
! --- parameters      
       
      TYPE(t_region), POINTER :: pRegion
      
! --- local variables
            
      INTEGER :: dLoc,ifcp,hLoc,nCells,vLoc      
      TYPE(t_global), POINTER :: global             
            
! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_EnforceSymmetryOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME   
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Enforcing optimal LES '// & 
                                 'integral symmetry...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer and variables
! ==============================================================================

      pGrid => pRegion%grid

      nCells = SIZE(pGrid%fsOLES,1)

! ==============================================================================
!     Enforce symmetry
! ==============================================================================          

! ------------------------------------------------------------------------------
!     Integral I4
! ------------------------------------------------------------------------------ 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 4...'
      END IF ! global%verbLevel  

      DO ifcp = 1,3
        DO vLoc = 1,9*nCells*nCells        
          IF ( symMapI45OLES(vLoc) /= 0 ) THEN 
            pGrid%int40OLES(:,ifcp,vLoc) = & 
              pGrid%int40OLES(:,ifcp,symMapI45OLES(vLoc))
            pGrid%int41OLES(:,ifcp,vLoc) = & 
              pGrid%int41OLES(:,ifcp,symMapI45OLES(vLoc))
            pGrid%int42OLES(:,ifcp,vLoc) = & 
              pGrid%int42OLES(:,ifcp,symMapI45OLES(vLoc))                
          END IF ! symMapI45OLES
        END DO ! vLoc
      END DO ! ifcp
   
! ------------------------------------------------------------------------------
!     Integral I5
! ------------------------------------------------------------------------------ 
      
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 5...'
      END IF ! global%verbLevel   
   
      DO ifcp = 1,3 ! equality of rows
        DO vLoc = 1,9*nCells*nCells
          IF ( symMapI45OLES(vLoc) /= 0 ) THEN 
            pGrid%int50OLES(ifcp,vLoc,:) = & 
              pGrid%int50OLES(ifcp,symMapI45OLES(vLoc),:)
            pGrid%int51OLES(ifcp,vLoc,:) = & 
              pGrid%int51OLES(ifcp,symMapI45OLES(vLoc),:)
            pGrid%int52OLES(ifcp,vLoc,:) = & 
              pGrid%int52OLES(ifcp,symMapI45OLES(vLoc),:)                
          END IF ! symMapI45OLES
        END DO ! vLoc
      END DO ! ifcp
      
! This is screwing up the matrix - no longer singular (yes, should be because of
! symmetry of weights). This is strange, but have not yet found out why destroying
! singularity. Enforcing the symmetries should make it as singular as possible...      
!      DO ifcp = 1,3 ! symmetry of entire matrix
!        DO dLoc = 1,9*nCells*nCells
!          DO vLoc = dLoc+1,9*nCells*nCells
!            hLoc = vLoc
!
!            pGrid%int50OLES(ifcp,dLoc,vLoc) = pGrid%int50OLES(ifcp,vLoc,dLoc)
!            pGrid%int51OLES(ifcp,dLoc,vLoc) = pGrid%int51OLES(ifcp,vLoc,dLoc)
!            pGrid%int52OLES(ifcp,dLoc,vLoc) = pGrid%int52OLES(ifcp,vLoc,dLoc)      
!          END DO ! hLoc            
!        END DO ! dLoc
!      END DO ! ifcp         
      
! ==============================================================================
!     End
! ==============================================================================
  
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Enforcing optimal LES '// & 
                                 'integral symmetry done.'
      END IF ! global%verbLevel   
  
      CALL DeregisterFunction(global)     
   
    END SUBROUTINE RFLU_EnforceSymmetryOLES    
    
    
    

! ******************************************************************************
!   Define the correlation functions which are to be integrated over the cells
! ******************************************************************************

! ==============================================================================
!   Second-order two-point correlation (for I^2), part 0
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation220(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: iFun,iFunNZ,j,l

! --- Start

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr22(iFunNZ)
        CALL RFLU_MapK2IJ(iFun,l,j)
        
        f(iFunNZ) = kd(l,j)
      END DO ! iFunNZ

! --- End 

    END SUBROUTINE RFLU_DefineCorrelation220

! ==============================================================================
!   Second-order two-point correlation (for I^2), part 1
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation221(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: iFun,iFunNZ,j,l
      DOUBLE PRECISION :: fact,rm2,rm23,rm2i

! --- Start

      rm2  = (z(4)-z(1))**2 + (z(5)-z(2))**2 + (z(6)-z(3))**2
      rm23 = rm2**(1.0D0/3.0D0)       
      rm2i = 1.0_RFREAL/MakeNonZero(rm2)

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr22(iFunNZ)
        CALL RFLU_MapK2IJ(iFun,l,j)
        
        f(iFunNZ) = rm23*(rm2i*(z(l+3)-z(l))*(z(j+3)-z(j)) - 4.0D0*kd(l,j))
      END DO ! iFunNZ

! --- End 

    END SUBROUTINE RFLU_DefineCorrelation221



! ==============================================================================
!   Third-order two-point correlation (for I^1)
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation32(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: i,iFun,iFunNZ,l,n
      DOUBLE PRECISION :: term1,term2,term3
      DOUBLE PRECISION :: zd(nDim+1)

! --- Start

      n = nzLoc

! --- Copy integration points and add constant face coordinate

      zd(1:nDim) = z(1:nDim)
      zd(nDim+1) = nzVal

! --- Evaluate integral

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr32(iFunNZ)
        CALL RFLU_MapK2IJ(iFun,l,i)
                         
        term1 = kd(i,n)*(zd(l)-zd(mapSurf2Vol2(n,l)))
        term2 = kd(i,l)*(zd(n)-zd(mapSurf2Vol2(n,n)))
        term3 = kd(n,l)*(zd(i)-zd(mapSurf2Vol2(n,i)))        
        
        f(iFunNZ) = nzSgn*(term1 - 1.5D0*(term2 + term3))
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation32



! ==============================================================================
!   Fourth-order three-point correlation (for I^4), part 0
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation430(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: i,iFun,iFunNZ,l,m,n
      DOUBLE PRECISION :: term1,term2,term3

! --- Start

      n = nzLoc

! --- Evaluate integral

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr43(iFunNZ)
        CALL RFLU_MapL2IJK(iFun,l,m,i)
        
        f(iFunNZ) = nzSgn*(kd(l,m)*kd(i,n) + kd(l,i)*kd(m,n) + kd(l,n)*kd(m,i))
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation430


! ==============================================================================
!   Fourth-order three-point correlation (for I^4), part 1
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation431(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: i,iFun,iFunNZ,j,k,l,m
      DOUBLE PRECISION :: r1m2,r1m23,r1m2i,r2m2,r2m23,r2m2i,r3m2,r3m23,r3m2i, & 
                          r4m2,r4m23,r4m2i,r5m2,r5m23,r5m2i,term1,term2,term3, & 
                          term4,term5
      DOUBLE PRECISION :: zd(nDim+1)

! --- Start

      j = nzLoc

! --- Copy integration points and add constant face coordinate

      zd(1:nDim) = z(1:nDim)
      zd(nDim+1) = nzVal

! --- Evaluate integral

      r1m2 = (zd(1) - zd(4))**2 & 
           + (zd(2) - zd(5))**2 & 
           + (zd(3) - zd(6))**2
      r2m2 = (zd(1) - zd(mapSurf2Vol3(j,1)))**2 & 
           + (zd(2) - zd(mapSurf2Vol3(j,2)))**2 & 
           + (zd(3) - zd(mapSurf2Vol3(j,3)))**2
      r3m2 = r2m2 
      r4m2 = (zd(4) - zd(mapSurf2Vol3(j,4)))**2 & 
           + (zd(5) - zd(mapSurf2Vol3(j,5)))**2 & 
           + (zd(6) - zd(mapSurf2Vol3(j,6)))**2           
      r5m2 = r4m2

      r1m23 = r1m2**(1.0D0/3.0D0)
      r2m23 = r2m2**(1.0D0/3.0D0)
      r3m23 = r2m23
      r4m23 = r4m2**(1.0D0/3.0D0)
      r5m23 = r4m23

      r1m2i = 1.0D0/MakeNonZero(r1m2)
      r2m2i = 1.0D0/MakeNonZero(r2m2)
      r3m2i = r2m2i
      r4m2i = 1.0D0/MakeNonZero(r4m2)
      r5m2i = r4m2i

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr43(iFunNZ)
        CALL RFLU_MapL2IJK(iFun,l,m,i)

        term1 = kd(i,j)*r1m23*(  r1m2i*(zd(l  )-zd(l+3)) & 
                                      *(zd(m  )-zd(m+3)) & 
                               - 4.0D0*kd(l,m))
        
        term2 = kd(l,i)*r5m23*(  r5m2i*(zd(m+3)-zd(mapSurf2Vol3(j,m))) & 
                                      *(zd(j+3)-zd(mapSurf2Vol3(j,j))) & 
                               - 4.0D0*kd(m,j))
                               
        term3 = kd(m,j)*r2m23*(  r2m2i*(zd(l  )-zd(mapSurf2Vol3(j,l))) & 
                                      *(zd(i  )-zd(mapSurf2Vol3(j,i))) & 
                               - 4.0D0*kd(l,i))
                               
        term4 = kd(l,j)*r4m23*(  r4m2i*(zd(m+3)-zd(mapSurf2Vol3(j,m))) & 
                                      *(zd(i+3)-zd(mapSurf2Vol3(j,i))) & 
                               - 4.0D0*kd(m,i))
                               
        term5 = kd(m,i)*r3m23*(  r3m2i*(zd(l  )-zd(mapSurf2Vol3(j,l))) & 
                                      *(zd(j  )-zd(mapSurf2Vol3(j,j))) & 
                               - 4.0D0*kd(l,j))
        f(iFunNZ) = nzSgn*(term1 + term2 + term3 + term4 + term5)        
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation431


! ==============================================================================
!   Fourth-order three-point correlation (for I^4), part 2
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation432(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: i,iFun,iFunNZ,j,k,l,m
      DOUBLE PRECISION :: r2m2,r2m23,r2m2i,r3m2,r3m23,r3m2i,r4m2,r4m23,r4m2i, & 
                          r5m2,r5m23,r5m2i,term1,term2
      DOUBLE PRECISION :: zd(nDim+1)

! --- Start

      j = nzLoc

! --- Copy integration points and add constant face coordinate

      zd(1:nDim) = z(1:nDim)
      zd(nDim+1) = nzVal

! --- Evaluate integral

      r2m2 = (zd(1) - zd(mapSurf2Vol3(j,1)))**2 & 
           + (zd(2) - zd(mapSurf2Vol3(j,2)))**2 & 
           + (zd(3) - zd(mapSurf2Vol3(j,3)))**2
      r3m2 = r2m2 
      r4m2 = (zd(4) - zd(mapSurf2Vol3(j,4)))**2 & 
           + (zd(5) - zd(mapSurf2Vol3(j,5)))**2 & 
           + (zd(6) - zd(mapSurf2Vol3(j,6)))**2           
      r5m2 = r4m2

      r2m23 = r2m2**(1.0D0/3.0D0)
      r3m23 = r2m23
      r4m23 = r4m2**(1.0D0/3.0D0)
      r5m23 = r4m23

      r2m2i = 1.0D0/MakeNonZero(r2m2)
      r3m2i = r2m2i
      r4m2i = 1.0D0/MakeNonZero(r4m2)
      r5m2i = r4m2i

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr43(iFunNZ)      
        CALL RFLU_MapL2IJK(iFun,l,m,i)
        
        term1 = & 
          r2m23*r5m23*(r2m2i*(zd(l  ) - zd(mapSurf2Vol3(j,l))) & 
                            *(zd(i  ) - zd(mapSurf2Vol3(j,i))) - 4.0D0*kd(l,i))             & 
                     *(r5m2i*(zd(m+3) - zd(mapSurf2Vol3(j,m))) & 
                            *(zd(j+3) - zd(mapSurf2Vol3(j,j))) - 4.0D0*kd(m,j))
        
        term2 = & 
          r3m23*r4m23*(r3m2i*(zd(l  ) - zd(mapSurf2Vol3(j,l))) & 
                            *(zd(j  ) - zd(mapSurf2Vol3(j,j))) - 4.0D0*kd(l,j))             & 
                     *(r4m2i*(zd(m+3) - zd(mapSurf2Vol3(j,m))) & 
                            *(zd(i+3) - zd(mapSurf2Vol3(j,i))) - 4.0D0*kd(m,i))
                               
        f(iFunNZ) = nzSgn*(term1 + term2)
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation432





! ==============================================================================
!   Fourth-order four-point correlation (for I^5), part 0
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation540(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: iFun,iFunNZ,j,l,k,m
      DOUBLE PRECISION :: term1,term2,term3

! --- Start, evaluate integral

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr44(iFunNZ)
        CALL RFLU_MapM2IJKL(iFun,l,m,j,k)
                
        f(iFunNZ) = kd(l,m)*kd(j,k) + kd(l,j)*kd(m,k) + kd(l,k)*kd(m,j)
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation540



! ==============================================================================
!   Fourth-order four-point correlation (for I^5), part 1
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation541(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: iFun,iFunNZ,j,k,l,m
      DOUBLE PRECISION :: r1m2,r1m23,r1m2i,r2m2,r2m23,r2m2i,r3m2,r3m23,r3m2i, & 
                          r4m2,r4m23,r4m2i,r5m2,r5m23,r5m2i,r6m2,r6m23,r6m2i, & 
                          term1,term2,term3,term4,term5,term6

! --- Start, evaluate integral

      r1m2 = (z(1) - z( 4))**2 & 
           + (z(2) - z( 5))**2 & 
           + (z(3) - z( 6))**2
      r2m2 = (z(1) - z( 7))**2 & 
           + (z(2) - z( 8))**2 & 
           + (z(3) - z( 9))**2
      r3m2 = (z(1) - z(10))**2 & 
           + (z(2) - z(11))**2 & 
           + (z(3) - z(12))**2
      r4m2 = (z(4) - z( 7))**2 & 
           + (z(5) - z( 8))**2 & 
           + (z(6) - z( 9))**2
      r5m2 = (z(4) - z(10))**2 & 
           + (z(5) - z(11))**2 & 
           + (z(6) - z(12))**2
      r6m2 = (z(7) - z(10))**2 & 
           + (z(8) - z(11))**2 & 
           + (z(9) - z(12))**2           

      r1m23 = r1m2**(1.0D0/3.0D0)
      r2m23 = r2m2**(1.0D0/3.0D0)
      r3m23 = r3m2**(1.0D0/3.0D0)
      r4m23 = r4m2**(1.0D0/3.0D0)
      r5m23 = r5m2**(1.0D0/3.0D0)
      r6m23 = r6m2**(1.0D0/3.0D0)

      r1m2i = 1.0D0/MakeNonZero(r1m2)
      r2m2i = 1.0D0/MakeNonZero(r2m2)
      r3m2i = 1.0D0/MakeNonZero(r3m2)
      r4m2i = 1.0D0/MakeNonZero(r4m2)
      r5m2i = 1.0D0/MakeNonZero(r5m2)
      r6m2i = 1.0D0/MakeNonZero(r6m2)

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr44(iFunNZ)              
        CALL RFLU_MapM2IJKL(iFun,l,m,j,k)
        
        term1 = r6m23*kd(l,m)*(r6m2i*(z(j+6) - z(j+9)) & 
                                    *(z(k+6) - z(k+9)) - 4.0D0*kd(j,k))
        term2 = r1m23*kd(j,k)*(r1m2i*(z(l  ) - z(l+3)) & 
                                    *(z(m  ) - z(m+3)) - 4.0D0*kd(l,m))
        term3 = r5m23*kd(l,j)*(r5m2i*(z(m+3) - z(m+9)) & 
                                    *(z(k+3) - z(k+9)) - 4.0D0*kd(m,k)) 
        term4 = r2m23*kd(m,k)*(r2m2i*(z(l  ) - z(l+6)) & 
                                    *(z(j  ) - z(j+6)) - 4.0D0*kd(l,j)) 
        term5 = r4m23*kd(l,k)*(r4m2i*(z(m+3) - z(m+6)) & 
                                    *(z(j+3) - z(j+6)) - 4.0D0*kd(m,j))                                     
        term6 = r3m23*kd(m,j)*(r3m2i*(z(l  ) - z(l+9)) & 
                                    *(z(k  ) - z(k+9)) - 4.0D0*kd(l,k))
                                                
        f(iFunNZ) = term1 + term2 + term3 + term4 + term5 + term6
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation541



! ==============================================================================
!   Fourth-order four-point correlation (for I^5), part 2
! ==============================================================================

    SUBROUTINE RFLU_DefineCorrelation542(nDim,z,nFunNZ,f)

! --- Arguments

      INTEGER, INTENT(IN) :: nDim,nFunNZ
      DOUBLE PRECISION, INTENT(IN) :: z(nDim)
      DOUBLE PRECISION, INTENT(INOUT) :: f(nFunNZ)

! --- Locals

      INTEGER :: iFun,iFunNZ,j,k,l,m
      DOUBLE PRECISION :: r1m2,r1m23,r1m2i,r2m2,r2m23,r2m2i,r3m2,r3m23,r3m2i, & 
                          r4m2,r4m23,r4m2i,r5m2,r5m23,r5m2i,r6m2,r6m23,r6m2i, & 
                          term1,term2,term3

! --- Start, evaluate integral

      r1m2 = (z(1) - z( 4))**2 & 
           + (z(2) - z( 5))**2 & 
           + (z(3) - z( 6))**2
      r2m2 = (z(1) - z( 7))**2 & 
           + (z(2) - z( 8))**2 & 
           + (z(3) - z( 9))**2
      r3m2 = (z(1) - z(10))**2 & 
           + (z(2) - z(11))**2 & 
           + (z(3) - z(12))**2
      r4m2 = (z(4) - z( 7))**2 & 
           + (z(5) - z( 8))**2 & 
           + (z(6) - z( 9))**2
      r5m2 = (z(4) - z(10))**2 & 
           + (z(5) - z(11))**2 & 
           + (z(6) - z(12))**2
      r6m2 = (z(7) - z(10))**2 & 
           + (z(8) - z(11))**2 & 
           + (z(9) - z(12))**2           

      r1m23 = r1m2**(1.0D0/3.0D0)
      r2m23 = r2m2**(1.0D0/3.0D0)
      r3m23 = r3m2**(1.0D0/3.0D0)
      r4m23 = r4m2**(1.0D0/3.0D0)
      r5m23 = r5m2**(1.0D0/3.0D0)
      r6m23 = r6m2**(1.0D0/3.0D0)

      r1m2i = 1.0D0/MakeNonZero(r1m2)
      r2m2i = 1.0D0/MakeNonZero(r2m2)
      r3m2i = 1.0D0/MakeNonZero(r3m2)
      r4m2i = 1.0D0/MakeNonZero(r4m2)
      r5m2i = 1.0D0/MakeNonZero(r5m2)
      r6m2i = 1.0D0/MakeNonZero(r6m2)

      DO iFunNZ = 1,nFunNZ
        iFun = mapFunNZ2FunCorr44(iFunNZ)              
        CALL RFLU_MapM2IJKL(iFun,l,m,j,k)
        
        term1 = & 
          r1m23*r6m23*(r1m2i*(z(l  ) - z(l+3))                  & 
                            *(z(m  ) - z(m+3)) - 4.0D0*kd(l,m)) & 
                     *(r6m2i*(z(j+6) - z(j+9))                  & 
                            *(z(k+6) - z(k+9)) - 4.0D0*kd(j,k))
        term2 = & 
          r2m23*r5m23*(r2m2i*(z(l  ) - z(l+6))                  & 
                            *(z(j  ) - z(j+6)) - 4.0D0*kd(l,j)) &
                     *(r5m2i*(z(m+3) - z(m+9))                  & 
                            *(z(k+3) - z(k+9)) - 4.0D0*kd(m,k))                             
        term3 = & 
          r3m23*r4m23*(r3m2i*(z(l  ) - z(l+9))                  & 
                            *(z(k  ) - z(k+9)) - 4.0D0*kd(l,k)) & 
                     *(r4m2i*(z(m+3) - z(m+6))                  & 
                            *(z(j+3) - z(j+6)) - 4.0D0*kd(m,j))          
                                                                        
        f(iFunNZ) = term1 + term2 + term3
      END DO ! iFunNZ      

! --- End

    END SUBROUTINE RFLU_DefineCorrelation542






! ******************************************************************************
!   Map from vector to second-order tensor
! ******************************************************************************

    SUBROUTINE RFLU_MapK2IJ(k,i,j)

! --- Arguments

      INTEGER, INTENT(IN) :: k
      INTEGER, INTENT(OUT) :: i,j

! --- Start

      j = (k+2)/3 ! NOTE integer division
      i = k - 3*(j-1)

! --- End

    END SUBROUTINE RFLU_MapK2IJ


! ******************************************************************************
!   Map from vector to third-order tensor
! ******************************************************************************

    SUBROUTINE RFLU_MapL2IJK(l,i,j,k)

! --- Arguments

      INTEGER, INTENT(IN) :: l
      INTEGER, INTENT(OUT) :: i,j,k

! --- Start

      k = (l+8)/9 ! NOTE integer division

      j = MOD((l+2)/3,3) ! NOTE integer division
      IF ( j == 0 ) THEN 
        j = 3
      END IF ! j
      
      i = l - 3*(j-1) - 9*(k-1)

! --- End

    END SUBROUTINE RFLU_MapL2IJK


! ******************************************************************************
!   Map from vector to fourth-order tensor
! ******************************************************************************

    SUBROUTINE RFLU_MapM2IJKL(m,i,j,k,l)

! --- Arguments

      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(OUT) :: i,j,k,l

! --- Start

      j = MOD((m+2)/3,3) ! NOTE integer division
      IF ( j == 0 ) THEN 
        j = 3
      END IF ! j
      
      k = MOD((m+8)/9,3) ! NOTE integer division
      IF ( k == 0 ) THEN 
        k = 3
      END IF ! k
      
      l = (m+26)/27 ! NOTE integer division           
      i = m - 27*(l-1) - 9*(k-1) - 3*(j-1)

! --- End

    END SUBROUTINE RFLU_MapM2IJKL



! ******************************************************************************
!   Get storage position for I1
! ******************************************************************************

    INTEGER FUNCTION RFLU_GetI1PosOLES(l,d)

! --- Arguments

      INTEGER, INTENT(IN) :: d,l

! --- Start

      RFLU_GetI1PosOLES = l + 3*(d-1)

! --- End

    END FUNCTION RFLU_GetI1PosOLES


! ******************************************************************************
!   Get storage position for I4
! ******************************************************************************

    INTEGER FUNCTION RFLU_GetI4PosOLES(l,m,d,e,nCells)

! --- Arguments

      INTEGER, INTENT(IN) :: d,e,l,m,nCells

! --- Start

      RFLU_GetI4PosOLES = m + 3*(l-1) + 9*(e-1) + 9*nCells*(d-1)

! --- End

    END FUNCTION RFLU_GetI4PosOLES



! ******************************************************************************
!   Get storage position for L
! ******************************************************************************

    INTEGER FUNCTION RFLU_GetLPosOLES(j,a)

! --- Arguments

      INTEGER, INTENT(IN) :: a,j

! --- Start

      RFLU_GetLPosOLES = j + 3*(a-1)

! --- End

    END FUNCTION RFLU_GetLPosOLES


! ******************************************************************************
!   Given storage position for L, find j and a 
! ******************************************************************************

    SUBROUTINE RFLU_GetLPosInvOLES(loc,j,a)

! --- Arguments

      INTEGER, INTENT(IN) :: loc
      INTEGER, INTENT(OUT) :: a,j

! --- Start

      a = (loc+2)/3
      j = loc - 3*(a-1)

! --- End

    END SUBROUTINE RFLU_GetLPosInvOLES
    


! ******************************************************************************
!   Get storage position for Q
! ******************************************************************************

    INTEGER FUNCTION RFLU_GetQPosOLES(j,k,b,g,nCells)

! --- Arguments

      INTEGER, INTENT(IN) :: b,g,j,k,nCells

! --- Start

      RFLU_GetQPosOLES = k + 3*(j-1) + 9*(g-1) + 9*nCells*(b-1)

! --- End

    END FUNCTION RFLU_GetQPosOLES
    
    
! ******************************************************************************
!   Given storage position for Q, find j,k, and b,g 
! ******************************************************************************

    SUBROUTINE RFLU_GetQPosInvOLES(loc,nCells,j,k,b,g)

! --- Arguments

      INTEGER, INTENT(IN) :: loc,nCells
      INTEGER, INTENT(OUT) :: b,g,j,k

! --- Start

      j = MOD((loc+2)/3,3)
      IF ( j == 0 ) THEN 
        j = 3 
      END IF ! j
  
      g = MOD((loc+8)/9,nCells) 
      IF ( g == 0 ) THEN 
        g = nCells
      END IF ! g

      b = (loc+9*nCells-1)/(9*nCells)            
            
      k = loc - 3*(j-1) - 9*(g-1) - 9*nCells*(b-1)

! --- End

    END SUBROUTINE RFLU_GetQPosInvOLES





! ******************************************************************************
!   Destroy geometry
! ******************************************************************************  

    SUBROUTINE RFLU_DestroyStencilsWeightsOLES(pRegion)
  
      IMPLICIT NONE
      
! --- parameters      
      
      TYPE(t_region), POINTER :: pRegion

! --- local variables

      INTEGER :: errorFlag
      TYPE(t_grid), POINTER :: pGrid              
      TYPE(t_global), POINTER :: global 

! ==============================================================================
!     Start
! ==============================================================================
      
      global => pRegion%global
      
      CALL RegisterFunction(global,'RFLU_DestroyStencilsWeightsOLES',&
  'RFLU_ModOLES.F90')

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying optimal LES '// & 
                                 'stencils and weights...'
      END IF ! global%verbLevel

! ==============================================================================
!     Set grid pointer
! ==============================================================================

      pGrid => pRegion%grid

! ==============================================================================
!     Deallocate memory
! ==============================================================================

      DEALLOCATE(pGrid%fsOLES,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%fsOLES')
      END IF ! global%error

      DEALLOCATE(pGrid%intLimOLES,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%intLimOLES')
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying optimal LES '// & 
                                 'stencils and weights done.'
      END IF ! global%verbLevel 
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_DestroyStencilsWeightsOLES




! ******************************************************************************
!   Allocate DCUHRE arrays
! ******************************************************************************

    SUBROUTINE RFLU_AllocateDCUHREArrays(global,nDim,nFunNZ,nFun)
    
      IMPLICIT NONE
      
! --- parameters      
      
      INTEGER, INTENT(IN) :: nDim,nFunNZ,nFun        
      TYPE(t_global), POINTER :: global 

! --- locals

      INTEGER :: errorFlag

! ==============================================================================
!     Start
! ==============================================================================
      
      CALL RegisterFunction(global,'RFLU_AllocateDCUHREArrays',&
  'RFLU_ModOLES.F90')

! ==============================================================================
!     Allocate memory
! ==============================================================================

      ALLOCATE(lowLim(nDim),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'lowLim')
      END IF ! global%error

      ALLOCATE(uppLim(nDim),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'uppLim')
      END IF ! global%error

      ALLOCATE(errAbsEst(nFun),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'errAbsEst')
      END IF ! global%error

      ALLOCATE(integralNZ(nFunNZ),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'integralNZ')
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      CALL DeregisterFunction(global)        
    
    END SUBROUTINE RFLU_AllocateDCUHREArrays





! ******************************************************************************
!   Deallocate DCUHRE arrays
! ******************************************************************************

    SUBROUTINE RFLU_DeallocateDCUHREArrays(global)
    
      IMPLICIT NONE
      
! --- arguments

      TYPE(t_global), POINTER :: global       
      
! --- locals
  
      INTEGER :: errorFlag      
      
! ==============================================================================
!     Start
! ==============================================================================
      
      CALL RegisterFunction(global,'RFLU_DeallocateDCUHREArrays',&
  'RFLU_ModOLES.F90')

! ==============================================================================
!     Allocate memory
! ==============================================================================

      DEALLOCATE(lowLim,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'lowLim')
      END IF ! global%error

      DEALLOCATE(uppLim,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'uppLim')
      END IF ! global%error

      DEALLOCATE(errAbsEst,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'errAbsEst')
      END IF ! global%error

      DEALLOCATE(integralNZ,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'integralNZ')
      END IF ! global%error

! ==============================================================================
!     End
! ==============================================================================

      CALL DeregisterFunction(global)        
    
    END SUBROUTINE RFLU_DeallocateDCUHREArrays




! ******************************************************************************
!   Set mapping for second-order two-point correlation (for I^2)
! ******************************************************************************

    SUBROUTINE RFLU_SetMapFunNZ2FunCorr22(nFunNZ)

! --- Arguments

      INTEGER, INTENT(IN) :: nFunNZ
      
! --- Locals

      INTEGER :: iFun,iFunNZ,j,l  
      INTEGER, PARAMETER :: NFUN = 9
      
! --- Start

      iFunNZ = 0

      IF ( nFunNZ /= NFUN ) THEN 
        DO iFun = 1,NFUN
          CALL RFLU_MapK2IJ(iFun,l,j)
          
          IF ( l == j ) THEN 
            iFunNZ = iFunNZ + 1
            mapFunNZ2FunCorr22(iFunNZ) = iFun
          END IF ! i
        END DO ! iFun
        
        mapFunNZ2FunCorr22(iFunNZ+1:NFUN) = CRAZY_VALUE_INT        
      ELSE 
        DO iFun = 1,NFUN
          iFunNZ = iFunNZ + 1
          mapFunNZ2FunCorr22(iFunNZ) = iFun
        END DO ! iFun
      END IF ! nFunNZ

! --- End
    
    END SUBROUTINE RFLU_SetMapFunNZ2FunCorr22    



! ******************************************************************************
!   Set mapping for third-order two-point correlation (for I^1)
! ******************************************************************************

    SUBROUTINE RFLU_SetMapFunNZ2FunCorr32(nFunNZ)

! --- Arguments

      INTEGER, INTENT(IN) :: nFunNZ
      
! --- Locals

      INTEGER :: i,iFun,iFunNZ,l  
      INTEGER, PARAMETER :: NFUN = 9
      
! --- Start

      iFunNZ = 0

      IF ( nFunNZ /= NFUN ) THEN 
        DO iFun = 1,NFUN
          CALL RFLU_MapK2IJ(iFun,l,i)
          
          IF ( l == i ) THEN 
            iFunNZ = iFunNZ + 1
            mapFunNZ2FunCorr32(iFunNZ) = iFun
          END IF ! i
        END DO ! iFun
        
        mapFunNZ2FunCorr32(iFunNZ+1:NFUN) = CRAZY_VALUE_INT        
      ELSE 
        DO iFun = 1,NFUN
          iFunNZ = iFunNZ + 1
          mapFunNZ2FunCorr32(iFunNZ) = iFun
        END DO ! iFun
      END IF ! nFunNZ

! --- End
    
    END SUBROUTINE RFLU_SetMapFunNZ2FunCorr32
    
   
    
    
! ******************************************************************************
!   Set mapping for fourth-order three-point correlation (for I^4)
! ******************************************************************************

    SUBROUTINE RFLU_SetMapFunNZ2FunCorr43(nFunNZ)

! --- Arguments

      INTEGER, INTENT(IN) :: nFunNZ
      
! --- Locals

      INTEGER :: i,iFun,iFunNZ,l,m,s  
      INTEGER, PARAMETER :: NFUN = 27
      REAL(RFREAL) :: term
      
! --- Start

      iFunNZ = 0
      
      s = nzLoc

      IF ( nFunNZ /= NFUN ) THEN 
        DO iFun = 1,NFUN
          CALL RFLU_MapL2IJK(iFun,l,m,i)
                    
          term = kd(l,m)*kd(i,s) + kd(l,i)*kd(m,s) + kd(l,s)*kd(m,i)
          
          IF ( NINT(term) /= 0 ) THEN 
            iFunNZ = iFunNZ + 1
            mapFunNZ2FunCorr43(iFunNZ) = iFun
          END IF ! i
        END DO ! iFun
        
        mapFunNZ2FunCorr43(iFunNZ+1:NFUN) = CRAZY_VALUE_INT        
      ELSE 
        DO iFun = 1,NFUN
          iFunNZ = iFunNZ + 1
          mapFunNZ2FunCorr43(iFunNZ) = iFun
        END DO ! iFun
      END IF ! nFunNZ

! --- End
    
    END SUBROUTINE RFLU_SetMapFunNZ2FunCorr43       




! ******************************************************************************
!   Set mapping for fourth-order four-point correlation (for I^5)
! ******************************************************************************

    SUBROUTINE RFLU_SetMapFunNZ2FunCorr44(nFunNZ)

! --- Arguments

      INTEGER, INTENT(IN) :: nFunNZ
      
! --- Locals

      INTEGER :: iFun,iFunNZ,j,k,l,m
      INTEGER, PARAMETER :: NFUN = 81
      REAL(RFREAL) :: term
      
! --- Start

      iFunNZ = 0
      
      IF ( nFunNZ /= NFUN ) THEN 
        DO iFun = 1,NFUN
          CALL RFLU_MapM2IJKL(iFun,l,m,j,k)
                    
          term = kd(l,m)*kd(j,k) + kd(l,j)*kd(m,k) + kd(l,k)*kd(m,j)
                    
          IF ( NINT(term) /= 0 ) THEN 
            iFunNZ = iFunNZ + 1
            mapFunNZ2FunCorr44(iFunNZ) = iFun
          END IF ! i
        END DO ! iFun
        
        mapFunNZ2FunCorr44(iFunNZ+1:NFUN) = CRAZY_VALUE_INT        
      ELSE 
        DO iFun = 1,NFUN
          iFunNZ = iFunNZ + 1
          mapFunNZ2FunCorr44(iFunNZ) = iFun
        END DO ! iFun
      END IF ! nFunNZ

! --- End
    
    END SUBROUTINE RFLU_SetMapFunNZ2FunCorr44       




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModOLES


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModOLES.F90,v $
!   Revision 1.10  2008/12/06 08:44:23  mtcampbe
!   Updated license.
!
!   Revision 1.9  2008/11/19 22:17:34  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.8  2004/01/22 16:03:59  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.7  2004/01/11 02:17:34  jiao
!   Eliminated some redundant trailing spaces that made some lines too long.
!   This changed was necessary to compile with NAG F90 compiler.
!
!   Revision 1.6  2003/05/16 02:27:44  haselbac
!   Removed KIND=RFREAL from NINT statements
!
!   Revision 1.5  2003/03/15 18:14:17  haselbac
!   Added KIND qualifyer to NINT statements
!
!   Revision 1.4  2002/10/08 15:49:21  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.3  2002/09/09 15:10:35  haselbac
!   global now under regions, many bug fixes as a result of CTR stay
!
!   Revision 1.2  2002/07/27 18:11:28  haselbac
!   Corrected I41, I42, I51, and I52 correlations, now get correct symmetries
!
!   Revision 1.1  2002/07/25 14:51:54  haselbac
!   Initial revision
!
! ******************************************************************************
















