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
! Purpose: Suite of routines to differentiate functions at boundary face 
!   centroids.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModDifferentiationBFaces.F90,v 1.6 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModDifferentiationBFaces

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI
  
  USE RFLU_ModConstraintUtils
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ComputeGradBFaces_1D, & 
            RFLU_ComputeGradBFaces, & 
            RFLU_ComputeGradBFacesWrapper, &
            RFLU_ComputeBFGradConstrWrapper
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModDifferentiationBFaces.F90,v $ $Revision: 1.6 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  


! ******************************************************************************
!
! Purpose: Compute 1D gradients of any vector or scalar at boundary face 
!   centers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeGradBFaces_1D(pRegion,pPatch,iBegVar,iEndVar, &
                                     iBegGrad,iEndGrad,var,grad)

  USE RFLU_ModPatchUtils, ONLY: RFLU_GetPatchNormalDirection
  USE RFLU_ModWeights, ONLY: RFLU_ComputeWtsX2C_1D

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: fnDirFlag,ifgIncludeFlag
  INTEGER :: errorFlag,fnDir,icg,ifg,iGrad,isl,iVar,nMembsMax,nMembs,order
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: locs,wts
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start, set pointers and variables
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeGradBFaces_1D',&
  'RFLU_ModDifferentiationBFaces.F90' )

#ifdef ROCPROF
  CALL FPROFILER_BEGINS("RFLU::ComputeGradBFaces_1D")
#endif

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! iEndVar

  pGrid => pRegion%grid

! ******************************************************************************    
! Compute gradients
! ******************************************************************************

  IF ( pPatch%bcType /= BC_VIRTUAL ) THEN
    nMembsMax = pPatch%bf2cs1DInfo%nCellMembsMax
    order     = 1 ! Order of derivative

    ifgIncludeFlag = .FALSE.

! ==============================================================================
!   Check patch geometry
! ==============================================================================

    IF ( pPatch%flatFlag .EQV. .FALSE. ) THEN 
      CALL ErrorStop(global,ERR_PATCH_NOT_FLAT,__LINE__)
    ELSE 
      CALL RFLU_GetPatchNormalDirection(global,pPatch,fnDir,fnDirFlag)

      IF ( fnDirFlag .EQV. .FALSE. ) THEN  
        CALL ErrorStop(global,ERR_PATCH_NOT_ALIGNED,__LINE__)
      END IF ! FloatEqual                          
    END IF ! pPatch%flatFlag    

! ==============================================================================
!   Allocate temporary memory
! ==============================================================================

    ALLOCATE(wts(0:nMembsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'wts')
    END IF ! global%error  

    ALLOCATE(locs(0:nMembsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'locs')
    END IF ! global%error   

! ==============================================================================
!   Loop over faces, compute gradients
! ==============================================================================

! ------------------------------------------------------------------------------
!   Include face ifg in stencil
! ------------------------------------------------------------------------------

    IF ( ifgIncludeFlag .EQV. .TRUE. ) THEN 
      DO ifg = 1,pPatch%nBFaces
        DO iGrad = iBegGrad,iEndGrad
          grad(XCOORD,iGrad,ifg) = 0.0_RFREAL
          grad(YCOORD,iGrad,ifg) = 0.0_RFREAL
          grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL                        
        END DO ! iGrad

        nMembs = pPatch%bf2cs1D(ifg)%nCellMembs

        locs(0) = pPatch%fc(fnDir,ifg)

        DO isl = 1,nMembs
          icg = pPatch%bf2cs1D(ifg)%cellMembs(isl)

          locs(isl) = pGrid%cofg(fnDir,icg)
        END DO ! isl

        CALL RFLU_ComputeWtsX2C_1D(global,order,nMembs+1,locs(0:nMembs), &
                                   pPatch%fc(fnDir,ifg),wts(0:nMembs))

        iGrad = iBegGrad

        DO iVar = iBegVar,iEndVar
! TEMPORARY
! Can only include ifg in stencil if have values on boundary. For the moment,
! that is not the case, but once NSCBC works, we can include values on
! boundary
!          grad(fnDir,iGrad,ifg) = wts(0)* 
! END TEMPORARY

          DO isl = 1,nMembs
            icg = pPatch%bf2cs1D(ifg)%cellMembs(isl)

            grad(fnDir,iGrad,ifg) = grad(fnDir,iGrad,ifg) &
                                  + wts(isl)*var(iVar,icg)
          END DO ! isl

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! ifg
           
! ------------------------------------------------------------------------------
!   Do not include face ifg in stencil
! ------------------------------------------------------------------------------
      
    ELSE
      DO ifg = 1,pPatch%nBFaces
        DO iGrad = iBegGrad,iEndGrad
          grad(XCOORD,iGrad,ifg) = 0.0_RFREAL
          grad(YCOORD,iGrad,ifg) = 0.0_RFREAL
          grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL                        
        END DO ! iGrad

        nMembs = pPatch%bf2cs1D(ifg)%nCellMembs

        DO isl = 1,nMembs
          icg = pPatch%bf2cs1D(ifg)%cellMembs(isl)

          locs(isl) = pGrid%cofg(fnDir,icg)
        END DO ! isl

        CALL RFLU_ComputeWtsX2C_1D(global,order,nMembs,locs(1:nMembs), &
                                   pPatch%fc(fnDir,ifg),wts(1:nMembs))

        iGrad = iBegGrad

        DO iVar = iBegVar,iEndVar
          DO isl = 1,nMembs
            icg = pPatch%bf2cs1D(ifg)%cellMembs(isl)

            grad(fnDir,iGrad,ifg) = grad(fnDir,iGrad,ifg) &
                                  + wts(isl)*var(iVar,icg)
          END DO ! isl

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! ifg  
    END IF ! ifgIncludeFlag

! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================

    DEALLOCATE(wts,STAT=errorFlag)                   
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'wts')
    END IF ! global%error
    
    DEALLOCATE(locs,STAT=errorFlag)                   
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'locs')
    END IF ! global%error    
  END IF ! pPatch%bcType

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF
  CALL FPROFILER_ENDS("RFLU::ComputeGradBFaces_1D")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeGradBFaces_1D






! ******************************************************************************
!
! Purpose: Compute gradients of any vector or scalar at boundary faces. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at boundary face centers
!
! Notes:
!   1. The face gradients differ from the cell gradients in that they are 
!      computed as weighted sums of variables rather than variable differences
!      because there are no variables located at faces.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeGradBFaces(pRegion,pPatch,iBegVar,iEndVar,iBegGrad, &
                                  iEndGrad,var,grad)

  USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_patch), POINTER :: pPatch  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,ifg,ifgBeg,ifgEnd,ifl,iGrad,isl,iVar
  REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz,r11, &
                  r12,r13,r14,r22,r23,r24,r33,r34,r44,term,term1,term2, &
                  term3,term4,wx,wy,wz
  REAL(RFREAL) :: fc(XCOORD:ZCOORD) 
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeGradBFaces',&
  'RFLU_ModDifferentiationBFaces.F90' )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! nVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Compute gradients
! ******************************************************************************

  IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
  
! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!       var(iVar,icg) =                                                  & 
!                     + REAL(4*(iVar-1)  ,RFREAL)                        & 
!                     - REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) &
!                     + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) &
!                     - REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG

! ==============================================================================   
!   Loop over faces and compute gradients 
!   Linear Interpolation Least Square Formulation
! ==============================================================================   

! ------------------------------------------------------------------------------
!   Select appropriate dimensionality
! ------------------------------------------------------------------------------

    SELECT CASE ( pRegion%mixtInput%dimens )

! --- Two dimensions ----------------------------------------------------------- 

      CASE ( 2 ) 

        DO ifl = 1,pPatch%nBFaces
          r11 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11)           
          r12 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_12) 
          r22 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_22)           
          r13 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_13) 
          r23 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_23)                   
          r33 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_33)                  

          c11 = 1.0_RFREAL/r11
          c22 = 1.0_RFREAL/r22
          c33 = 1.0_RFREAL/r33

          c12 = - c11*r12
          c13 = -(c11*r13 + c12*c22*r23)  

          c23 = - c22*r23

          fc(XCOORD) = pPatch%fc(XCOORD,ifl)
          fc(YCOORD) = pPatch%fc(YCOORD,ifl) 

          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,ifl) = 0.0_RFREAL
            grad(YCOORD,iGrad,ifl) = 0.0_RFREAL
            grad(ZCOORD,iGrad,ifl) = 0.0_RFREAL
          END DO ! iVar

          DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
            icg = pPatch%bf2cs(ifl)%cellMembs(isl)

            dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
            dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)   

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

            dx = term*dx
            dy = term*dy        

            term1 = c11*c11*(                    dx)
            term2 = c22*c22*(           dy + c12*dx)
            term3 = c33*c33*(term + c23*dy + c13*dx)                           

! TEMPORARY 
            wx = term*(term1 + c12*term2 + c13*term3)
            wy = term*(            term2 + c23*term3)

!            wx = term*(term1 + c12*term2)
!            wy = term*(            term2)
! END TEMPORARY

            iGrad = iBegGrad

            DO iVar = iBegVar,iEndVar              
! TEMPORARY 
              grad(XCOORD,iGrad,ifl) = grad(XCOORD,iGrad,ifl) + wx*var(iVar,icg)
              grad(YCOORD,iGrad,ifl) = grad(YCOORD,iGrad,ifl) + wy*var(iVar,icg)
 
!              dVar = var(iVar,icg) - pPatch%mixt%cv(iVar,ifl)	           
!              grad(XCOORD,iGrad,ifl) = grad(XCOORD,iGrad,ifl) + wx*dVar
!              grad(YCOORD,iGrad,ifl) = grad(YCOORD,iGrad,ifl) + wy*dVar
! END TEMPORARY

              iGrad = iGrad + 1
            END DO ! iVar
          END DO ! isl        

        END DO ! ifl  

! --- Three dimensions ---------------------------------------------------------

      CASE ( 3 ) 
        DO ifl = 1,pPatch%nBFaces                

          r11 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11)           
          r12 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_12) 
          r22 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_22)           
          r13 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_13) 
          r23 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_23)                   
          r33 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_33)
          r14 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_14)
          r24 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_24)
          r34 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_34)
          r44 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_44)                  

          c11 = 1.0_RFREAL/r11
          c22 = 1.0_RFREAL/r22
          c33 = 1.0_RFREAL/r33
          c44 = 1.0_RFREAL/r44

          c12 = - c11*r12
          c13 = -(c11*r13 + c12*c22*r23) 
          c14 = -(c11*r14 + c12*c22*r24 + c13*c33*r34) 

          c23 = - c22*r23
          c24 = -(c22*r24 + c23*c33*r34)

          c34 = - c33*r34

          fc(XCOORD) = pPatch%fc(XCOORD,ifl)
          fc(YCOORD) = pPatch%fc(YCOORD,ifl) 
          fc(ZCOORD) = pPatch%fc(ZCOORD,ifl)             

          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,ifl) = 0.0_RFREAL
            grad(YCOORD,iGrad,ifl) = 0.0_RFREAL
            grad(ZCOORD,iGrad,ifl) = 0.0_RFREAL
          END DO ! iVar

          DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
            icg = pPatch%bf2cs(ifl)%cellMembs(isl)

            dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
            dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
            dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)   

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

            dx = term*dx
            dy = term*dy
            dz = term*dz        

            term1 = c11*c11*(                             dx)
            term2 = c22*c22*(                    dy + c12*dx)
            term3 = c33*c33*(           dz + c23*dy + c13*dx)
            term4 = c44*c44*(term + c34*dz + c24*dy + c14*dx)               

            wx = term*(term1 + c12*term2 + c13*term3 + c14*term4)
            wy = term*(            term2 + c23*term3 + c24*term4)
            wz = term*(                        term3 + c34*term4)

            iGrad = iBegGrad

            DO iVar = iBegVar,iEndVar              
              grad(XCOORD,iGrad,ifl) = grad(XCOORD,iGrad,ifl) + wx*var(iVar,icg)
              grad(YCOORD,iGrad,ifl) = grad(YCOORD,iGrad,ifl) + wy*var(iVar,icg)
              grad(ZCOORD,iGrad,ifl) = grad(ZCOORD,iGrad,ifl) + wz*var(iVar,icg)

              iGrad = iGrad + 1
            END DO ! iVar              
          END DO ! isl

        END DO ! ifl

! --- Default ------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)        
    END SELECT ! pRegion%mixtInput%dimens

! DEBUG
!    ifgBeg = 1
!    ifgEnd = pPatch%nBFaces
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MINVAL(grad(YCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MINVAL(grad(ZCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MAXVAL(grad(XCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MAXVAL(grad(YCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MAXVAL(grad(ZCOORD,iGrad,ifgBeg:ifgEnd))
!    END DO ! iGrad
! END DEBUG

  END IF ! pPatch%bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeGradBFaces









! ******************************************************************************
!
! Purpose: Compute constrained gradients of any vector or scalar at boundary 
!   faces. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   varInfo     Variable information
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at boundary face centers
!
! Notes:
!   1. The face gradients differ from the cell gradients in that they are 
!      computed as weighted sums of variables rather than variable differences
!      because there are no variables located at faces.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeBFGradConstr(pRegion,pPatch,iBegVar,iEndVar,iBegGrad, & 
                                    iEndGrad,varInfo,var,grad)

  USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  INTEGER, INTENT(IN) :: varInfo(iBegVar:iEndVar)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_patch), POINTER :: pPatch  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,iCol,ifg,ifgBeg,ifgEnd,ifl,ifl2,iGrad, & 
             iPatch2,iRow,isl,iVar,ix,iy,iz,nCols,nConstr,nRows,sCount
  INTEGER, DIMENSION(:), ALLOCATABLE :: constrType
  REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz,gx,gy, & 
                  gz,r11,r12,r13,r14,r22,r23,r24,r33,r34,r44,term,term1, &
                  term2,term3,term4,varf,varc    
  REAL(RFREAL) :: colMax(4)
  REAL(RFREAL) :: fc(XCOORD:ZCOORD)    
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch2

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeBFGradConstr',&
  'RFLU_ModDifferentiationBFaces.F90' )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! nVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Compute gradients
! ******************************************************************************

  IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
  
! DEBUG
!  DO icg = 1,pGrid%nCellsTot
!    DO iVar = iBegVar,iEndVar
!      var(iVar,icg) =                                                  & 
!                    + REAL(4*(iVar-1)  ,RFREAL)                        & 
!                    - REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) &
!                    + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) 
!                    - REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!    END DO ! iVar
!  END DO ! icg
! END DEBUG

! ==============================================================================   
!   Loop over faces and compute gradients 
! ==============================================================================   

    DO ifl = 1,pPatch%nBFaces

! ------------------------------------------------------------------------------
!     Initialize gradients
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,ifl) = 0.0_RFREAL
        grad(YCOORD,iGrad,ifl) = 0.0_RFREAL
        grad(ZCOORD,iGrad,ifl) = 0.0_RFREAL
      END DO ! iGrad

      fc(XCOORD) = pPatch%fc(XCOORD,ifl)
      fc(YCOORD) = pPatch%fc(YCOORD,ifl)
      fc(ZCOORD) = pPatch%fc(ZCOORD,ifl)            

      ALLOCATE(constrType(0:pPatch%bf2cs(ifl)%nBFaceMembs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'constrType')
      END IF ! global%error                       

! ------------------------------------------------------------------------------
!     Compute gradients 
! ------------------------------------------------------------------------------

      iGrad = iBegGrad

      DO iVar = iBegVar,iEndVar

! ----- Determine whether bface itself is constrained -------------------------- 

        constrType(0) = RFLU_GetConstrType(pRegion,pPatch,varInfo(iVar),ifl)          

        IF ( constrType(0) /= CONSTR_TYPE_DIRICHLET ) THEN 
          constrType(0) = CONSTR_TYPE_NONE
          
          varf = 0.0_RFREAL ! IMPORTANT
        ELSE
          varf = RFLU_GetConstrValue(pRegion,pPatch,varInfo(iVar),ifl) 
        END IF ! constrType                  

! ----- Determine number of constraints ----------------------------------------

        nConstr = 0       

        DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs
          iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
          ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

          pPatch2 => pRegion%patches(iPatch2)

          constrType(isl) = RFLU_GetConstrType(pRegion,pPatch2,varInfo(iVar), &
                                               ifl2)

          IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
            nConstr = nConstr + 1
          ELSE 
            constrType(isl) = CONSTR_TYPE_NONE
          END IF ! constrType          
        END DO ! isl          

! ------------------------------------------------------------------------------
!       Gradients constrained by Dirichlet boundary conditions. Treated as 
!       soft constraints. NOTE do not need to treat case of unconstrained
!       gradients here because they should already have been computed. If 
!       they have not been computed at this stage, they will be set to zero
!       by initialization above.
! ------------------------------------------------------------------------------      

        IF ( (constrType(0) == CONSTR_TYPE_DIRICHLET) .OR. (nConstr > 0) ) THEN 

! ------- Allocate temporary memory 

          nRows = pPatch%bf2cs(ifl)%nCellMembs + nConstr

          SELECT CASE ( constrType(0) )
            CASE ( CONSTR_TYPE_DIRICHLET )
              nCols = pRegion%mixtInput%dimens
            CASE ( CONSTR_TYPE_NONE ) 
              nCols = pRegion%mixtInput%dimens + 1
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! constrType(0)

          ALLOCATE(a(nRows,nCols),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'a')
          END IF ! global%error

          ALLOCATE(aInv(nCols,nRows),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'aInv')
          END IF ! global%error

! ------- Define left-hand side matrix 

          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 2 )    
              SELECT CASE ( constrType(0) )
                CASE ( CONSTR_TYPE_DIRICHLET )
                  DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
                    icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                    dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                    dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)

                    term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                    a(isl,1) = term*dx
                    a(isl,2) = term*dy
                  END DO ! isl  

                  iRow = pPatch%bf2cs(ifl)%nCellMembs              

                  DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs
                    IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                      iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                      ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                      pPatch2 => pRegion%patches(iPatch2)

                      dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                      dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)

                      term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                      iRow = iRow + 1

                      a(iRow,1) = term*dx
                      a(iRow,2) = term*dy
                    END IF ! constrType
                  END DO ! isl                              
                CASE ( CONSTR_TYPE_NONE ) 
                  DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
                    icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                    dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                    dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)

                    term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                    a(isl,1) = term
                    a(isl,2) = term*dx
                    a(isl,3) = term*dy
                  END DO ! isl  

                  iRow = pPatch%bf2cs(ifl)%nCellMembs              

                  DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs
                    IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                      iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                      ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                      pPatch2 => pRegion%patches(iPatch2)

                      dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                      dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)

                      term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                      iRow = iRow + 1

                      a(iRow,1) = term
                      a(iRow,2) = term*dx
                      a(iRow,3) = term*dy
                    END IF ! constrType
                  END DO ! isl                      
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! constrType(0)                                                                                                       
            CASE ( 3 )
              SELECT CASE ( constrType(0) )
                CASE ( CONSTR_TYPE_DIRICHLET )
                  DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
                    icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                    dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                    dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
                    dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)               

                    term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                    a(isl,1) = term*dx
                    a(isl,2) = term*dy
                    a(isl,3) = term*dz
                  END DO ! isl  

                  iRow = pPatch%bf2cs(ifl)%nCellMembs 

                  DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs              
                    IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                      iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                      ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                      pPatch2 => pRegion%patches(iPatch2)

                      dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                      dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)
                      dz = pPatch2%fc(ZCOORD,ifl2) - fc(ZCOORD)                  

                      term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                      iRow = iRow + 1

                      a(iRow,1) = term*dx
                      a(iRow,2) = term*dy
                      a(iRow,3) = term*dz
                    END IF ! constrType                 
                  END DO ! isl                                                              
                CASE ( CONSTR_TYPE_NONE ) 
                  DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
                    icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                    dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                    dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
                    dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)               

                    term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                    a(isl,1) = term
                    a(isl,2) = term*dx
                    a(isl,3) = term*dy
                    a(isl,4) = term*dz
                  END DO ! isl  

                  iRow = pPatch%bf2cs(ifl)%nCellMembs 

                  DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs              
                    IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                      iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                      ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                      pPatch2 => pRegion%patches(iPatch2)

                      dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                      dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)
                      dz = pPatch2%fc(ZCOORD,ifl2) - fc(ZCOORD)                  

                      term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                      iRow = iRow + 1

                      a(iRow,1) = term
                      a(iRow,2) = term*dx
                      a(iRow,3) = term*dy
                      a(iRow,4) = term*dz
                    END IF ! constrType                 
                  END DO ! isl                                      
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! constrType(0)                                     
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
          END SELECT ! pRegion%mixtInput%dimens 

! ------- Compute constrained gradient weights 

          DO iCol = 1,nCols          
            colMax(iCol) = -HUGE(1.0_RFREAL)

            DO iRow = 1,nRows
              colMax(iCol) = MAX(colMax(iCol),ABS(a(iRow,iCol)))
            END DO ! iRow

            DO iRow = 1,nRows
              a(iRow,iCol) = a(iRow,iCol)/colMax(iCol)
            END DO ! iRow                     
          END DO ! iCol                          

          CALL RFLU_InvertMatrixSVD(global,nRows,nCols,a,aInv,sCount)

          DO iCol = 1,nCols
            DO iRow = 1,nRows
              aInv(iCol,iRow) = aInv(iCol,iRow)/colMax(iCol) 
            END DO ! iRow
          END DO ! iCol                         

! TEMPORARY
          IF ( sCount /= 0 ) THEN
            WRITE(*,*) 'ERROR - Singular matrix in RFLU_ComputeGradBFacesConstr!'
            STOP
          END IF ! sCount
! END TEMPORARY

! ------- Compute constrained gradients 

          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 2 )
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL

              SELECT CASE ( constrType(0) )
                CASE ( CONSTR_TYPE_DIRICHLET )
                  ix = 1
                  iy = 2            
                CASE ( CONSTR_TYPE_NONE ) 
                  ix = 2
                  iy = 3
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
              END SELECT ! constrType(0)
              
              DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs
                icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy) 

                gx = gx + term*aInv(ix,isl)*(var(iVar,icg) - varf)
                gy = gy + term*aInv(iy,isl)*(var(iVar,icg) - varf)
              END DO ! isl                    

              iRow = pPatch%bf2cs(ifl)%nCellMembs

              DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN                   
                  iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                  ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                  pPatch2 => pRegion%patches(iPatch2)

                  dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)

                  term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                  varc = RFLU_GetConstrValue(pRegion,pPatch2,varInfo(iVar),ifl2)                

                  iRow = iRow + 1

                  gx = gx + term*aInv(ix,iRow)*(varc - varf)
                  gy = gy + term*aInv(iy,iRow)*(varc - varf)
                END IF ! constrType     
              END DO ! isl

              grad(XCOORD,iGrad,ifl) = gx
              grad(YCOORD,iGrad,ifl) = gy
              grad(ZCOORD,iGrad,ifl) = 0.0_RFREAL           

            CASE ( 3 ) 
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL
              gz = 0.0_RFREAL    

              SELECT CASE ( constrType(0) )
                CASE ( CONSTR_TYPE_DIRICHLET )
                  ix = 1
                  iy = 2 
                  iz = 3           
                CASE ( CONSTR_TYPE_NONE ) 
                  ix = 2
                  iy = 3
                  iz = 4
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
              END SELECT ! constrType(0)

              DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs
                icg = pPatch%bf2cs(ifl)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
                dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)               

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                gx = gx + term*aInv(ix,isl)*(var(iVar,icg) - varf)
                gy = gy + term*aInv(iy,isl)*(var(iVar,icg) - varf)
                gz = gz + term*aInv(iz,isl)*(var(iVar,icg) - varf)              
              END DO ! isl                    

              iRow = pPatch%bf2cs(ifl)%nCellMembs

              DO isl = 1,pPatch%bf2cs(ifl)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN               
                  iPatch2 = pPatch%bf2cs(ifl)%bFaceMembs(1,isl)        
                  ifl2    = pPatch%bf2cs(ifl)%bFaceMembs(2,isl)        

                  pPatch2 => pRegion%patches(iPatch2)

                  dx = pPatch2%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch2%fc(YCOORD,ifl2) - fc(YCOORD)
                  dz = pPatch2%fc(ZCOORD,ifl2) - fc(ZCOORD)               

                  term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)   

                  varc = RFLU_GetConstrValue(pRegion,pPatch2,varInfo(iVar),ifl2)

                  iRow = iRow + 1

                  gx = gx + term*aInv(ix,iRow)*(varc - varf)
                  gy = gy + term*aInv(iy,iRow)*(varc - varf)
                  gz = gz + term*aInv(iz,iRow)*(varc - varf)
                END IF ! constrType                        
              END DO ! isl

              grad(XCOORD,iGrad,ifl) = gx
              grad(YCOORD,iGrad,ifl) = gy
              grad(ZCOORD,iGrad,ifl) = gz                                            
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
          END SELECT ! pRegion%mixtInput%dimens   

! ------- Deallocate temporary memory          

          DEALLOCATE(a,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'a')
          END IF ! global%error

          DEALLOCATE(aInv,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'aInv')
          END IF ! global%error            
        END IF ! nConstr

        iGrad = iGrad + 1                                  
      END DO ! iVar      

      DEALLOCATE(constrType,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'constrType')
      END IF ! global%error              

    END DO ! ifl

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,ifgBeg:ifgEnd)), & 
!                       MINVAL(grad(YCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MINVAL(grad(ZCOORD,iGrad,ifgBeg:ifgEnd)), & 
!                       MAXVAL(grad(XCOORD,iGrad,ifgBeg:ifgEnd)), & 
!                       MAXVAL(grad(YCOORD,iGrad,ifgBeg:ifgEnd)), &
!                       MAXVAL(grad(ZCOORD,iGrad,ifgBeg:ifgEnd))
!    END DO ! iGrad
! END DEBUG

  END IF ! pPatch%bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeBFGradConstr


  
  
  
  
! ******************************************************************************
!
! Purpose: Compute gradients of any vector or scalar at boundary face centroids. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeBFGradConstrWrapper(pRegion,pPatch,iBegVar,iEndVar, & 
                                           iBegGrad,iEndGrad,varInfo,var,grad)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  INTEGER, INTENT(IN) :: varInfo(iBegVar:iEndVar)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_patch), POINTER :: pPatch  
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeBFGradConstrWrapper',&
  'RFLU_ModDifferentiationBFaces.F90' )

! ******************************************************************************    
! Call gradient routines
! ******************************************************************************    

  SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
    CASE ( 1 ) 
! TO DO
! END TO DO         
    CASE ( 2,3 ) 
      CALL RFLU_ComputeBFGradConstr(pRegion,pPatch,iBegVar,iEndVar, & 
                                    iBegGrad,iEndGrad,varInfo,var,grad)                                         
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pMixtInput%stencilDimensBFaces

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeBFGradConstrWrapper  
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Compute gradients of any vector or scalar at boundary face centroids. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradBFacesWrapper(pRegion,pPatch,iBegVar,iEndVar, &
                                           iBegGrad,iEndGrad,var,grad)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradBFacesWrapper',&
  'RFLU_ModDifferentiationBFaces.F90' )

! ******************************************************************************    
!   Call gradient routines
! ******************************************************************************    

    SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
      CASE ( 1 )
        CALL RFLU_ComputeGradBFaces_1D(pRegion,pPatch,iBegVar,iEndVar, &
                                       iBegGrad,iEndGrad,var,grad)               
      CASE ( 2,3 ) 
        CALL RFLU_ComputeGradBFaces(pRegion,pPatch,iBegVar,iEndVar,iBegGrad, &
                                    iEndGrad,var,grad)                                         
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradBFacesWrapper
  
  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModDifferentiationBFaces


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDifferentiationBFaces.F90,v $
! Revision 1.6  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:39:05  mparmar
! Removed bf2bg
!
! Revision 1.3  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.2  2006/04/07 14:46:21  haselbac
! Rewrite in terms of wrapper funcs bcos of 1D routines
!
! Revision 1.1  2005/10/27 19:31:35  haselbac
! Initial revision
!
! ******************************************************************************
  











