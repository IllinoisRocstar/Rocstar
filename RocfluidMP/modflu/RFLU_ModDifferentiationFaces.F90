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
! Purpose: Suite of routines to differentiate functions at face centroids.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModDifferentiationFaces.F90,v 1.9 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModDifferentiationFaces

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
  PUBLIC :: RFLU_ComputeGradFacesWrapper, &
            RFLU_ComputeGradFacesConstr, &
            RFLU_ComputeGradConstrained
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModDifferentiationFaces.F90,v $ $Revision: 1.9 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  





! ******************************************************************************
!
! Purpose: Compute gradients of any vector or scalar at face centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at face centers
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
    
  SUBROUTINE RFLU_ComputeGradFaces(pRegion,iBegVar,iEndVar,iBegGrad,iEndGrad, &
                                   var,grad)

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
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,ifg,ifl,iGrad,iPatch,isl,iVar
    REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz,r11, &
                    r12,r13,r14,r22,r23,r24,r33,r34,r44,term,term1,term2, & 
                    term3,term4,wx,wy,wz
    REAL(RFREAL) :: fc(XCOORD:ZCOORD)    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradFaces',&
  'RFLU_ModDifferentiationFaces.F90')

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradFaces")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************    
!   Loop over faces and compute gradients 
! ******************************************************************************    

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!       var(iVar,icg) =                                                  &
!                     + REAL(4*(iVar-1)  ,RFREAL)                        & 
!                     + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) &
!                     + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) &  
!                     + REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!     END DO ! iVar
!   END DO ! icg
! END DEBUG
                  
! ==============================================================================
!   Select appropriate dimensionality
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%dimens )
        
! ------------------------------------------------------------------------------
!     Two dimensions
! ------------------------------------------------------------------------------
         
      CASE ( 2 ) 
        DO ifg = 1,pGrid%nFaces    
          r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
          r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
          r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
          r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
          r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
          r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)                  

          c11 = 1.0_RFREAL/r11
          c22 = 1.0_RFREAL/r22
          c33 = 1.0_RFREAL/r33

          c12 = - c11*r12
          c13 = -(c11*r13 + c12*c22*r23) 

          c23 = - c22*r23

          fc(XCOORD) = pGrid%fc(XCOORD,ifg)
          fc(YCOORD) = pGrid%fc(YCOORD,ifg)                    

          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,ifg) = 0.0_RFREAL
            grad(YCOORD,iGrad,ifg) = 0.0_RFREAL
            grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL
          END DO ! iVar

          DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
            icg = pGrid%f2cs(ifg)%cellMembs(isl)

            dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
            dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

            dx = term*dx
            dy = term*dy       

            term1 = c11*c11*(                    dx)
            term2 = c22*c22*(           dy + c12*dx)
            term3 = c33*c33*(term + c23*dy + c13*dx)                                 

            wx = term*(term1 + c12*term2 + c13*term3)
            wy = term*(            term2 + c23*term3)

            iGrad = iBegGrad
            
            DO iVar = iBegVar,iEndVar
              grad(XCOORD,iGrad,ifg) = grad(XCOORD,iGrad,ifg) + wx*var(iVar,icg)
              grad(YCOORD,iGrad,ifg) = grad(YCOORD,iGrad,ifg) + wy*var(iVar,icg)

              iGrad = iGrad + 1          
            END DO ! iVar                        
          END DO ! isl
        END DO ! ifg

! ------------------------------------------------------------------------------
!      Three dimensions
! ------------------------------------------------------------------------------

      CASE ( 3 ) 
        DO ifg = 1,pGrid%nFaces        
          r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
          r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
          r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
          r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
          r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
          r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)
          r14 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_14)
          r24 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_24)
          r34 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_34)
          r44 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_44)                  

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

          fc(XCOORD) = pGrid%fc(XCOORD,ifg)
          fc(YCOORD) = pGrid%fc(YCOORD,ifg) 
          fc(ZCOORD) = pGrid%fc(ZCOORD,ifg)                               

          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,ifg) = 0.0_RFREAL
            grad(YCOORD,iGrad,ifg) = 0.0_RFREAL
            grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL
          END DO ! iVar

          DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
            icg = pGrid%f2cs(ifg)%cellMembs(isl)

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
              grad(XCOORD,iGrad,ifg) = grad(XCOORD,iGrad,ifg) + wx*var(iVar,icg)
              grad(YCOORD,iGrad,ifg) = grad(YCOORD,iGrad,ifg) + wy*var(iVar,icg)
              grad(ZCOORD,iGrad,ifg) = grad(ZCOORD,iGrad,ifg) + wz*var(iVar,icg)

              iGrad = iGrad + 1          
            END DO ! iVar
          END DO ! isl
        END DO ! ifg

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------
        
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MINVAL(grad(YCOORD,iGrad,1:pGrid%nFaces)), &
!                       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nFaces)), &
!                       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nFaces))
!    END DO ! iGrad    
! END DEBUG

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradFaces")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradFaces








! ******************************************************************************
!
! Purpose: Compute constrained gradients of any vector or scalar at face 
!   centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   varInfo     Variable information
!   var         Variables of which gradients are to be determined
!
! Output:
!   grad        Gradients of variables at face centers
!
! Notes:
!   1. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   2. Restricted to linear reconstruction for now.
!   3. Restricted to Dirichlet constraints for now.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradFacesConstr(pRegion,iBegVar,iEndVar,iBegGrad, &
                                         iEndGrad,varInfo,var,grad)

    USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
    INTEGER, INTENT(IN) :: varInfo(iBegVar:iEndVar)
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,iCol,ifg,ifl,ifl2,iGrad,iPatch,isl,iRow,iVar, & 
               nCols,nConstr,nRows,sCount
    INTEGER, DIMENSION(:), ALLOCATABLE :: constrType               
    REAL(RFREAL) :: cwt,dx,dy,dz,term,varc,wtx,wty,wtz,gx,gy,gz
    REAL(RFREAL) :: colMax(4)
    REAL(RFREAL) :: fc(XCOORD:ZCOORD) 
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradFacesConstr',&
  'RFLU_ModDifferentiationFaces.F90')

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradFacesConstr")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    cwt = pRegion%mixtInput%cReconstFacesWeight

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!       var(iVar,icg) =                                                  &
!                     + REAL(4*(iVar-1)  ,RFREAL)                        &  
!                     + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) &
!                     + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) &  
!                     + REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG

! ******************************************************************************    
!   Loop over faces and compute gradients 
! ******************************************************************************    

    DO ifl = 1,pGrid%nFacesConstr
      ifg = pGrid%ifgConstr(ifl)

! ==============================================================================
!     Initialize gradients
! ==============================================================================

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,ifg) = 0.0_RFREAL
        grad(YCOORD,iGrad,ifg) = 0.0_RFREAL
        grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL
      END DO ! iGrad
      
      fc(XCOORD) = pGrid%fc(XCOORD,ifg)
      fc(YCOORD) = pGrid%fc(YCOORD,ifg)
      fc(ZCOORD) = pGrid%fc(ZCOORD,ifg)            
     
      ALLOCATE(constrType(pGrid%f2cs(ifg)%nBFaceMembs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'constrType')
      END IF ! global%error       
      
! ==============================================================================
!     Compute gradients 
! ==============================================================================

      iGrad = iBegGrad

      DO iVar = iBegVar,iEndVar
        
! ------------------------------------------------------------------------------
!       Determine number of constraints
! ------------------------------------------------------------------------------  
              
        nConstr = 0       
              
        DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
          iPatch = pGrid%f2cs(ifg)%bFaceMembs(1,isl)        
          ifl2   = pGrid%f2cs(ifg)%bFaceMembs(2,isl)        
         
          pPatch => pRegion%patches(iPatch)
         
          constrType(isl) = RFLU_GetConstrType(pRegion,pPatch,varInfo(iVar), &
                                               ifl2)
                          
          IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
            nConstr = nConstr + 1
          ELSE 
            constrType(isl) = CONSTR_TYPE_NONE
          END IF ! constrType          
        END DO ! pGrid%f2cs(ifg)%nBFaceMembs 

! ------------------------------------------------------------------------------
!       Gradients constrained by Dirichlet boundary conditions. Treated as 
!       soft constraints. NOTE do not need to treat case of unconstrained
!       gradients here because they should already have been computed. If 
!       they have not been computed at this stage, they will be set to zero
!       by initialization above.
! ------------------------------------------------------------------------------      
                 
        IF ( nConstr > 0 ) THEN 

! ------- Allocate temporary memory --------------------------------------------

          nRows = pGrid%f2cs(ifg)%nCellMembs + nConstr
          nCols = pRegion%mixtInput%dimens + 1

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

! ------- Define left-hand side matrix -----------------------------------------
 
          SELECT CASE ( pRegion%mixtInput%dimens ) 

! --------- Two dimensions          

            CASE ( 2 )             
              DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                a(isl,1) = term
                a(isl,2) = term*dx
                a(isl,3) = term*dy
              END DO ! isl  
                            
              iRow = pGrid%f2cs(ifg)%nCellMembs
              
              DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                  iPatch = pGrid%f2cs(ifg)%bFaceMembs(1,isl)        
                  ifl2   = pGrid%f2cs(ifg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl2) - fc(YCOORD)

                  term = cwt/SQRT(dx*dx + dy*dy)            

                  iRow = iRow + 1
                  
                  a(iRow,1) = term
                  a(iRow,2) = term*dx
                  a(iRow,3) = term*dy
                END IF ! constrType
              END DO ! isl              

! --------- Three dimensions              
                            
            CASE ( 3 )
              DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
                dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)                

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                a(isl,1) = term
                a(isl,2) = term*dx
                a(isl,3) = term*dy
                a(isl,4) = term*dz
              END DO ! isl  
              
              iRow = pGrid%f2cs(ifg)%nCellMembs 
              
              DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                  iPatch = pGrid%f2cs(ifg)%bFaceMembs(1,isl)        
                  ifl2   = pGrid%f2cs(ifg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl2) - fc(YCOORD)
                  dz = pPatch%fc(ZCOORD,ifl2) - fc(ZCOORD)                  

                  term = cwt/SQRT(dx*dx + dy*dy + dz*dz)            

                  iRow = iRow + 1

                  a(iRow,1) = term
                  a(iRow,2) = term*dx
                  a(iRow,3) = term*dy
                  a(iRow,4) = term*dz
                END IF ! constrType                 
              END DO ! isl                  
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
          END SELECT ! pRegion%mixtInput%dimens 

! ------- Compute constrained gradient weights ---------------------------------
                                 
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
            WRITE(*,*) 'ERROR - Singular matrix in RFLU_ComputeGradFacesConstr!'
            STOP
          END IF ! sCount     
! END TEMPORARY
   
! ------- Compute constrained gradient weights ---------------------------------

          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 2 )
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL
    
              DO isl = 1,pGrid%f2cs(ifg)%nCellMembs
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
 
                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy) 
 
                gx = gx + term*aInv(2,isl)*var(iVar,icg)
                gy = gy + term*aInv(3,isl)*var(iVar,icg)
              END DO ! isl                    
              
              iRow = pGrid%f2cs(ifg)%nCellMembs

              DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN                   
                  iPatch = pGrid%f2cs(ifg)%bFaceMembs(1,isl)        
                  ifl2   = pGrid%f2cs(ifg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl2) - fc(YCOORD)

                  term = cwt/SQRT(dx*dx + dy*dy)            

                  varc = RFLU_GetConstrValue(pRegion,pPatch,varInfo(iVar),ifl2)                

                  iRow = iRow + 1

                  gx = gx + term*aInv(2,iRow)*varc
                  gy = gy + term*aInv(3,iRow)*varc
                END IF ! constrType     
              END DO ! isl
                            
              grad(XCOORD,iGrad,ifg) = gx
              grad(YCOORD,iGrad,ifg) = gy
              grad(ZCOORD,iGrad,ifg) = 0.0_RFREAL                           
            CASE ( 3 ) 
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL
              gz = 0.0_RFREAL
    
              DO isl = 1,pGrid%f2cs(ifg)%nCellMembs
                icg = pGrid%f2cs(ifg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg) - fc(XCOORD)
                dy = pGrid%cofg(YCOORD,icg) - fc(YCOORD)
                dz = pGrid%cofg(ZCOORD,icg) - fc(ZCOORD)                
 
                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                gx = gx + term*aInv(2,isl)*var(iVar,icg)
                gy = gy + term*aInv(3,isl)*var(iVar,icg)
                gz = gz + term*aInv(4,isl)*var(iVar,icg)               
              END DO ! isl                    
              
              iRow = pGrid%f2cs(ifg)%nCellMembs

              DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN               
                  iPatch = pGrid%f2cs(ifg)%bFaceMembs(1,isl)        
                  ifl2   = pGrid%f2cs(ifg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl2) - fc(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl2) - fc(YCOORD)
                  dz = pPatch%fc(ZCOORD,ifl2) - fc(ZCOORD)                

                  term = cwt/SQRT(dx*dx + dy*dy + dz*dz)   

                  varc = RFLU_GetConstrValue(pRegion,pPatch,varInfo(iVar),ifl2)

                  iRow = iRow + 1

                  gx = gx + term*aInv(2,iRow)*varc
                  gy = gy + term*aInv(3,iRow)*varc
                  gz = gz + term*aInv(4,iRow)*varc
                END IF ! constrType                        
              END DO ! isl
              
              grad(XCOORD,iGrad,ifg) = gx
              grad(YCOORD,iGrad,ifg) = gy
              grad(ZCOORD,iGrad,ifg) = gz                                               
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
          END SELECT ! pRegion%mixtInput%dimens   
  
! ------- Deallocate temporary memory ------------------------------------------         

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
        END IF ! gradType

        iGrad = iGrad + 1                                  
      END DO ! iVar      
 
      DEALLOCATE(constrType,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'constrType')
      END IF ! global%error  
    END DO ! ifl

! ******************************************************************************
!   End
! ******************************************************************************

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MINVAL(grad(YCOORD,iGrad,1:pGrid%nFaces)), &
!                       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nFaces)), & 
!                       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nFaces)), &
!                       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nFaces))
!    END DO ! iGrad   
! END DEBUG

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradFacesConstr")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradFacesConstr







! ******************************************************************************
!
! Purpose: Wrapper routine for computing gradients of any vector or scalar at 
!   face centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
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
    
  SUBROUTINE RFLU_ComputeGradFacesWrapper(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,var,grad)

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

    CALL RegisterFunction(global,'RFLU_ComputeGradFacesWrapper',&
  'RFLU_ModDifferentiationFaces.F90' )

! ******************************************************************************    
!   Call gradient routines
! ******************************************************************************    

    SELECT CASE ( pRegion%mixtInput%stencilDimensFaces )
      CASE ( 1 ) 
! TO DO 
! END TO DO                       
      CASE ( 2,3 ) 
        CALL RFLU_ComputeGradFaces(pRegion,iBegVar,iEndVar,iBegGrad,iEndGrad, &
                                   var,grad)                                          
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradFacesWrapper






! ******************************************************************************
!
! Purpose: Compute constrained gradients. 
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   dimens              Dimensionality
!   nCellMembs          Number of cell members in stencil
!   nBFaceMembs         Number of boundary face members in stencil
!   order               Order of gradient reconstruction
!   dra                 Coordinate differences of cell members
!   drb                 Coordinate differences of boundary face members
!   rhsa                Variable values of cell members
!   rhsb                Variable values of boundary face members
!
! Output:
!   gradLocal           Gradients of variables
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeGradConstrained(global,dimens,nCellMembs,nBFaceMembs, &
                                         order,dra,drb,rhsa,rhsb,gradLocal)

    USE RFLU_ModStencilsUtils, ONLY: RFLU_ComputeStencilSize
    USE ModTools, ONLY: CompFact
    
    IMPLICIT NONE 
  
! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: dimens
    INTEGER :: nBFaceMembs,nCellMembs,order    
    REAL(RFREAL) :: gradLocal(XCOORD:XYZMAG),rhsa(nCellMembs),rhsb(nBFaceMembs)
    REAL(RFREAL) :: dra(XCOORD:ZCOORD,nCellMembs),drb(XCOORD:ZCOORD,nBFaceMembs)  
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,isl,j,lda,ldb,nCols,p,q,r,term,workArrayRealSize
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: workArrayReal
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,b

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeGradConstrained',&
  'RFLU_ModDifferentiationFaces.F90')

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nCols = RFLU_ComputeStencilSize(global,dimens,1,order)
    
    workArrayRealSize = nBFaceMembs + nCols + 100*nCellMembs

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(a(nCellMembs,nCols),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'a')
    END IF ! global%error     

    ALLOCATE(b(nBFaceMembs,nCols),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'b')
    END IF ! global%error

    ALLOCATE(workArrayReal(workArrayRealSize),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'workArrayReal')
    END IF ! global%error    

! ******************************************************************************
!   Build matrices
! ******************************************************************************

! ==============================================================================
!   Matrix arising from unconstrained part
! ==============================================================================

    DO isl = 1,nCellMembs
      a(isl,nCols) = 1.0_RFREAL

      j = 1

      DO p = 1,order
        DO q = 0,p
          DO r = 0,q
            term = 1.0_RFREAL/(CompFact(p-q)*CompFact(q-r)*CompFact(r))                 
            a(isl,j) = term*dra(XCOORD,isl)**(p-q) &
                           *dra(YCOORD,isl)**(q-r) &
                           *dra(ZCOORD,isl)**r 
            j = j + 1
          END DO ! r
        END DO ! q
      END DO ! p
    END DO ! isl    

! ==============================================================================
!   Matrix arising from constraints
! ==============================================================================

    DO isl = 1,nBFaceMembs
      b(isl,nCols) = 1.0_RFREAL

      j = 1

      DO p = 1,order
        DO q = 0,p
          DO r = 0,q
            term = 1.0_RFREAL/(CompFact(p-q)*CompFact(q-r)*CompFact(r))                 
            b(isl,j) = term*drb(XCOORD,isl)**(p-q) &
                           *drb(YCOORD,isl)**(q-r) &
                           *drb(ZCOORD,isl)**r 
            j = j + 1
          END DO ! r
        END DO ! q
      END DO ! p
    END DO ! isl  

! ******************************************************************************
!   Compute gradients
! ******************************************************************************

    lda = MAX(1,nCellMembs)
    ldb = MAX(1,nBFaceMembs)

    CALL dgglse(nCellMembs,nCols,nBFaceMembs,a,lda,b,ldb,rhsa,rhsb,gradLocal, & 
                workArrayReal,workArrayRealSize,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_LAPACK_OUTPUT,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(a,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'a')
    END IF ! global%error     

    DEALLOCATE(b,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'b')
    END IF ! global%error

    DEALLOCATE(workArrayReal,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'workArrayReal')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_ComputeGradConstrained






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModDifferentiationFaces


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDifferentiationFaces.F90,v $
! Revision 1.9  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.6  2006/04/07 14:47:03  haselbac
! Added wrapper func in anticipation of 1D routines
!
! Revision 1.5  2005/12/25 15:32:10  haselbac
! Added user-specified constraint weight
!
! Revision 1.4  2005/10/27 19:12:38  haselbac
! Bug fixes, separated constr and unconstr grads, bgrads now in own module
!
! Revision 1.3  2005/10/18 03:00:26  haselbac
! Bug fix: Incorrect computation of boundary gradients
!
! Revision 1.2  2005/10/07 20:02:23  haselbac
! Bug fix: Incorrect loop limits, lead to core dump for parallel cases
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************










