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
! Purpose: Suite of routines to differentiate function values at cell centroids.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModDifferentiationCells.F90,v 1.12 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModDifferentiationCells

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
  PUBLIC :: RFLU_ComputeGradCellsWrapper
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModDifferentiationCells.F90,v $ $Revision: 1.12 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  



! ******************************************************************************
!
! Purpose: Compute 1D gradients of any vector or scalar at cell centers.
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
    
  SUBROUTINE RFLU_ComputeGradCells_1D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                      iEndGrad,var,grad)

    USE RFLU_ModWeights, ONLY: RFLU_ComputeWtsX2C_1D

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

    LOGICAL :: icgIncludeFlag
    INTEGER :: errorFlag,fnDir,fnDirEnd,icg,icg2,iGrad,isl,iVar,nMembsMax, &
               nMembs,order
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: locs,wts
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCells_1D',&
  'RFLU_ModDifferentiationCells.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCells_1D")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    nMembsMax = pGrid%c2csInfo%nCellMembsMax
    order     = 1 ! Order of derivative

    icgIncludeFlag = .TRUE.

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        fnDirEnd = XCOORD
      CASE ( 2 )
        fnDirEnd = YCOORD
      CASE ( 3 )       
        fnDirEnd = ZCOORD        
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************    
!   Allocate temporary memory
! ******************************************************************************

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

! ******************************************************************************
!   Compute gradients
! ******************************************************************************

! ==============================================================================
!   Include cell icg in stencil
! ==============================================================================

    IF ( icgIncludeFlag .EQV. .TRUE. ) THEN 
      DO icg = 1,pGrid%nCellsTot
        DO iGrad = iBegGrad,iEndGrad
          grad(XCOORD,iGrad,icg) = 0.0_RFREAL
          grad(YCOORD,iGrad,icg) = 0.0_RFREAL
          grad(ZCOORD,iGrad,icg) = 0.0_RFREAL                        
        END DO ! iGrad

        DO fnDir = XCOORD,fnDirEnd
          nMembs = pGrid%c2cs1D(fnDir,icg)%nCellMembs

          locs(0) = pGrid%cofg(fnDir,icg)

          DO isl = 1,nMembs
            icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)

            locs(isl) = pGrid%cofg(fnDir,icg2)
          END DO ! isl

          CALL RFLU_ComputeWtsX2C_1D(global,order,nMembs+1,locs(0:nMembs), &
                                     pGrid%cofg(fnDir,icg),wts(0:nMembs))

          iGrad = iBegGrad

          DO iVar = iBegVar,iEndVar
            grad(fnDir,iGrad,icg) = wts(0)*var(iVar,icg)

            DO isl = 1,nMembs
              icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)

              grad(fnDir,iGrad,icg) = grad(fnDir,iGrad,icg) &
                                    + wts(isl)*var(iVar,icg2)
            END DO ! isl

            iGrad = iGrad + 1          
          END DO ! iVar
        END DO ! fnDir
      END DO ! icg

! ==============================================================================
!   Do not include cell icg in stencil
! ==============================================================================

    ELSE 
      DO icg = 1,pGrid%nCellsTot
        DO iGrad = iBegGrad,iEndGrad
          grad(XCOORD,iGrad,icg) = 0.0_RFREAL
          grad(YCOORD,iGrad,icg) = 0.0_RFREAL
          grad(ZCOORD,iGrad,icg) = 0.0_RFREAL                        
        END DO ! iGrad

        DO fnDir = XCOORD,fnDirEnd
          nMembs = pGrid%c2cs1D(fnDir,icg)%nCellMembs

          DO isl = 1,nMembs
            icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)

            locs(isl) = pGrid%cofg(fnDir,icg2)
          END DO ! isl

          CALL RFLU_ComputeWtsX2C_1D(global,order,nMembs,locs(1:nMembs), &
                                     pGrid%cofg(fnDir,icg),wts(1:nMembs))

          iGrad = iBegGrad

          DO iVar = iBegVar,iEndVar
            DO isl = 1,nMembs
              icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)

              grad(fnDir,iGrad,icg) = grad(fnDir,iGrad,icg) &
                                    + wts(isl)*var(iVar,icg2)
            END DO ! isl

            iGrad = iGrad + 1          
          END DO ! iVar
        END DO ! fnDir
      END DO ! icg      
    END IF ! icgIncludeFlag      
    
! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

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

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCells_1D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCells_1D







! ******************************************************************************
!
! Purpose: Compute 2D gradients of any vector or scalar at cell centers. 
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
! Notes: 
!   1. Use original form of reconstruction, i.e., weights multiplying variable
!      differences, because central location is cell, so have variable located
!      there. This is in contrast to face gradients and vertex interpolation,
!      where need to use modified form of reconstruction.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   3. Restricted to linear reconstruction for now.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCells_2D(pRegion,iBegVar,iEndVar,iBegGrad, &
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

    INTEGER :: errorFlag,icg,icg2,iGrad,isl,iVar
    REAL(RFREAL) :: c11,c12,c22,dVar,dx,dy,r11,r12,r22,term,term1,term2,wx,wy
    REAL(RFREAL) :: cofg(XCOORD:YCOORD)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCells_2D',&
  'RFLU_ModDifferentiationCells.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCells_2D")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************    
!   Loop over cells and compute gradients 
! ******************************************************************************    

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!	var(iVar,icg) = REAL(4*(iVar-1)  ,RFREAL)			 & 
!		      + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) & 
!		      + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG

    DO icg = 1,pGrid%nCellsTot                               
      r11 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11)           
      r12 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_12) 
      r22 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_22)           

      c11 = 1.0_RFREAL/r11
      c22 = 1.0_RFREAL/r22

      c12 = -c11*r12

      cofg(XCOORD) = pGrid%cofg(XCOORD,icg)
      cofg(YCOORD) = pGrid%cofg(YCOORD,icg)                    

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = 0.0_RFREAL
        grad(YCOORD,iGrad,icg) = 0.0_RFREAL
        grad(ZCOORD,iGrad,icg) = 0.0_RFREAL ! Set z-component to zero
      END DO ! iVar

      DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
        icg2 = pGrid%c2cs(icg)%cellMembs(isl)

        dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
        dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

        dx = term*dx
        dy = term*dy        

        term1 = c11*c11*(         dx)
        term2 = c22*c22*(dy + c12*dx)

        wx = term*(term1 + c12*term2)
        wy = term*(            term2)

        iGrad = iBegGrad       

        DO iVar = iBegVar,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wy*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl                   
    END DO ! icg 

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		        MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		        MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		        MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCells_2D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCells_2D







! ******************************************************************************
!
! Purpose: Compute 3D gradients of any vector or scalar at cell centers. 
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
! Notes: 
!   1. Use original form of reconstruction, i.e., weights multiplying variable
!      differences, because central location is cell, so have variable located
!      there. This is in contrast to face gradients and vertex interpolation,
!      where need to use modified form of reconstruction.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   3. Restricted to linear reconstruction for now.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCells_3D(pRegion,iBegVar,iEndVar,iBegGrad, &
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

    INTEGER :: errorFlag,icg,icg2,iGrad,isl,iVar
    REAL(RFREAL) :: c11,c12,c13,c22,c23,c33,dVar,dx,dy,dz,r11,r12,r13,r22, &
                    r23,r33,term,term1,term2,term3,wx,wy,wz
    REAL(RFREAL) :: cofg(XCOORD:ZCOORD)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCells_3D',&
  'RFLU_ModDifferentiationCells.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCells_3D")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************    
!   Loop over cells and compute gradients 
! ******************************************************************************    

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!        var(iVar,icg) = REAL(4*(iVar-1)  ,RFREAL)                        & 
!                      + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) & 
!                      + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) &
!                      + REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG
    
    DO icg = 1,pGrid%nCellsTot                       
      r11 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11)           
      r12 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_12) 
      r22 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_22)           
      r13 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_13) 
      r23 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_23)                   
      r33 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_33)

      c11 = 1.0_RFREAL/r11
      c22 = 1.0_RFREAL/r22
      c33 = 1.0_RFREAL/r33

      c12 = - c11*r12
      c13 = -(c11*r13 + c12*c22*r23) 

      c23 = - c22*r23

      cofg(XCOORD) = pGrid%cofg(XCOORD,icg)
      cofg(YCOORD) = pGrid%cofg(YCOORD,icg)
      cofg(ZCOORD) = pGrid%cofg(ZCOORD,icg)            

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = 0.0_RFREAL
        grad(YCOORD,iGrad,icg) = 0.0_RFREAL
        grad(ZCOORD,iGrad,icg) = 0.0_RFREAL
      END DO ! iVar

      DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
        icg2 = pGrid%c2cs(icg)%cellMembs(isl)

        dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
        dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)
        dz = pGrid%cofg(ZCOORD,icg2) - cofg(ZCOORD)   

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

        dx = term*dx
        dy = term*dy
        dz = term*dz        

        term1 = c11*c11*(                  dx)
        term2 = c22*c22*(         dy + c12*dx)
        term3 = c33*c33*(dz + c23*dy + c13*dx)

        wx = term*(term1 + c12*term2 + c13*term3)
        wy = term*(            term2 + c23*term3)
        wz = term*(                        term3)

        iGrad = iBegGrad       

        DO iVar = iBegVar,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wy*dVar
          grad(ZCOORD,iGrad,icg) = grad(ZCOORD,iGrad,icg) + wz*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl                   
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!                       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!                       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCells_3D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCells_3D








! ******************************************************************************
!
! Purpose: Compute 2D gradients of any vector or scalar at cell centers. 
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
! Notes: 
!   1. Use original form of reconstruction, i.e., weights multiplying variable
!      differences, because central location is cell, so have variable located
!      there. This is in contrast to face gradients and vertex interpolation,
!      where need to use modified form of reconstruction.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   3. Optimized by Adam Moody and Charles Shereda, LLNL.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCellsFast_2D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,var,grad)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegVar,iEndVar,iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: endisl,errorFlag,icg,icg2,icg3,icg4,iGrad,isl,iVar
    REAL(RFREAL) :: c11,c1111,c12,c22,c2222,dVar,dVar3,dVar4,dx,dx_icg,dy, &
                    dy_icg,r11,r12,r22,term,term1,term2,wtx,wty
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCellsFast_2D',&
  'RFLU_ModDifferentiationCells.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCellsFast_2D")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************    
!   Loop over cells and compute gradients 
! ******************************************************************************    

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!	var(iVar,icg) = REAL(4*(iVar-1)  ,RFREAL)			 & 
!		      + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) & 
!		      + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG

    DO icg = 1,pGrid%nCellsTot                       
      r11 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11)           
      r12 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_12) 
      r22 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_22)           

      c11 = 1.0_RFREAL/r11
      c22 = 1.0_RFREAL/r22

      c12 = -c11*r12

      endisl = pGrid%c2cs(icg)%nCellMembs 

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = 0.0_RFREAL
        grad(YCOORD,iGrad,icg) = 0.0_RFREAL
	grad(ZCOORD,iGrad,icg) = 0.0_RFREAL ! Set z-component to zero
      END DO ! iGrad

      icg3 = pGrid%c2cs(icg)%cellMembs(1)
      icg4 = pGrid%c2cs(icg)%cellMembs(2)
      
      dx_icg = pGrid%cofg(XCOORD,icg)
      dy_icg = pGrid%cofg(YCOORD,icg)

      dVar3 = var(iBegVar,icg3) - var(iBegVar,icg)
      dVar4 = var(iBegVar,icg4) - var(iBegVar,icg)

      c1111 = c11*c11
      c2222 = c22*c22

      DO isl = 1,endisl-2
        icg2 = icg3
        icg3 = icg4
        icg4 = pGrid%c2cs(icg)%cellMembs(isl+2)

        dVar  = dVar3
        dVar3 = dVar4
        dVar4 = var(iBegVar,icg4) - var(iBegVar,icg)

        dx = pGrid%cofg(XCOORD,icg2) - dx_icg
        dy = pGrid%cofg(YCOORD,icg2) - dy_icg

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

        dx = term*dx
        dy = term*dy

        term1 = c1111*(         dx)
        term2 = c2222*(dy + c12*dx)

        wtx = term*(term1 + c12*term2)
        wty = term*(            term2)

        iGrad = iBegGrad       

        grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
        grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar

        iGrad = iGrad + 1          
        
        DO iVar = iBegVar+1,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl    
                     
      DO isl = endisl-1,endisl
        icg2 = icg3
        icg3 = icg4

        dVar  = dVar3
        dVar3 = dVar4

        dx = pGrid%cofg(XCOORD,icg2) - dx_icg
        dy = pGrid%cofg(YCOORD,icg2) - dy_icg

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

        dx = term*dx
        dy = term*dy

        term1 = c1111*(         dx)
        term2 = c2222*(dy + c12*dx)

        wtx = term*(term1 + c12*term2)
        wty = term*(            term2)

        iGrad = iBegGrad       

        grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
        grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar

        iGrad = iGrad + 1          
        
        DO iVar = iBegVar+1,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCellsFast_2D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCellsFast_2D







! ******************************************************************************
!
! Purpose: Compute 3D gradients of any vector or scalar at cell centers. 
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
! Notes: 
!   1. Use original form of reconstruction, i.e., weights multiplying variable
!      differences, because central location is cell, so have variable located
!      there. This is in contrast to face gradients and vertex interpolation,
!      where need to use modified form of reconstruction.
!   2. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   3. Optimized by Adam Moody and Charles Shereda, LLNL.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCellsFast_3D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,var,grad)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegVar,iEndVar,iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:), POINTER :: var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: endisl,errorFlag,icg,icg2,icg3,icg4,iGrad,isl,iVar
    REAL(RFREAL) :: c11,c1111,c12,c13,c22,c2222,c23,c33,c3333,dVar,dVar3, &
                    dVar4,dx,dx_icg,dy,dy_icg,dz,dz_icg,r11,r12,r13,r22, &
                    r23,r33,term,term1,term2,term3,wtx,wty,wtz
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCellsFast_3D',&
  'RFLU_ModDifferentiationCells.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCellsFast_3D")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************    
!   Loop over cells and compute gradients 
! ******************************************************************************    

! DEBUG
!    DO icg = 1,pGrid%nCellsTot
!      DO iVar = iBegVar,iEndVar
!        var(iVar,icg) = REAL(4*(iVar-1)  ,RFREAL)                        & 
!                      + REAL(4*(iVar-1)+1,RFREAL)*pGrid%cofg(XCOORD,icg) & 
!                      + REAL(4*(iVar-1)+2,RFREAL)*pGrid%cofg(YCOORD,icg) &
!                      + REAL(4*(iVar-1)+3,RFREAL)*pGrid%cofg(ZCOORD,icg) 
!      END DO ! iVar
!    END DO ! icg
! END DEBUG

    DO icg = 1,pGrid%nCellsTot                       
      r11 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_11)           
      r12 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_12) 
      r22 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_22)           
      r13 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_13) 
      r23 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_23)                   
      r33 = pGrid%c2cs(icg)%xyzMoms(XYZ_MOM_33)

      c11 = 1.0_RFREAL/r11
      c22 = 1.0_RFREAL/r22
      c33 = 1.0_RFREAL/r33

      c12 = - c11*r12
      c13 = -(c11*r13 + c12*c22*r23) 

      c23 = - c22*r23

      endisl = pGrid%c2cs(icg)%nCellMembs 

      DO iGrad = iBegGrad,iEndGrad
	grad(XCOORD,iGrad,icg) = 0.0_RFREAL
	grad(YCOORD,iGrad,icg) = 0.0_RFREAL
	grad(ZCOORD,iGrad,icg) = 0.0_RFREAL
      END DO ! iGrad

      icg3 = pGrid%c2cs(icg)%cellMembs(1)
      icg4 = pGrid%c2cs(icg)%cellMembs(2)
      
      dx_icg = pGrid%cofg(XCOORD,icg)
      dy_icg = pGrid%cofg(YCOORD,icg)
      dz_icg = pGrid%cofg(ZCOORD,icg)

      dVar3 = var(iBegVar,icg3) - var(iBegVar,icg)
      dVar4 = var(iBegVar,icg4) - var(iBegVar,icg)

      c1111 = c11*c11
      c2222 = c22*c22
      c3333 = c33*c33

      DO isl = 1,endisl-2
        icg2 = icg3
        icg3 = icg4
        icg4 = pGrid%c2cs(icg)%cellMembs(isl+2)

        dVar  = dVar3
        dVar3 = dVar4
        dVar4 = var(iBegVar,icg4) - var(iBegVar,icg)

        dx = pGrid%cofg(XCOORD,icg2) - dx_icg
        dy = pGrid%cofg(YCOORD,icg2) - dy_icg
        dz = pGrid%cofg(ZCOORD,icg2) - dz_icg

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

        dx = term*dx
        dy = term*dy
        dz = term*dz        

        term1 = c1111*(                  dx)
        term2 = c2222*(         dy + c12*dx)
        term3 = c3333*(dz + c23*dy + c13*dx)

        wtx = term*(term1 + c12*term2 + c13*term3)
        wty = term*(            term2 + c23*term3)
        wtz = term*(                        term3)

        iGrad = iBegGrad       

        grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
        grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar
        grad(ZCOORD,iGrad,icg) = grad(ZCOORD,iGrad,icg) + wtz*dVar
	
        iGrad = iGrad + 1          
        
        DO iVar = iBegVar+1,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar
          grad(ZCOORD,iGrad,icg) = grad(ZCOORD,iGrad,icg) + wtz*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl    
                     
      DO isl = endisl-1,endisl
        icg2 = icg3
        icg3 = icg4

        dVar  = dVar3
        dVar3 = dVar4

        dx = pGrid%cofg(XCOORD,icg2) - dx_icg
        dy = pGrid%cofg(YCOORD,icg2) - dy_icg
        dz = pGrid%cofg(ZCOORD,icg2) - dz_icg

        term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

        dx = term*dx
        dy = term*dy
        dz = term*dz        

        term1 = c1111*(                  dx)
        term2 = c2222*(         dy + c12*dx)
        term3 = c3333*(dz + c23*dy + c13*dx)

        wtx = term*(term1 + c12*term2 + c13*term3)
        wty = term*(            term2 + c23*term3)
        wtz = term*(                        term3)

        iGrad = iBegGrad       

        grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
        grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar
        grad(ZCOORD,iGrad,icg) = grad(ZCOORD,iGrad,icg) + wtz*dVar

        iGrad = iGrad + 1          
        
        DO iVar = iBegVar+1,iEndVar
          dVar = var(iVar,icg2) - var(iVar,icg)

          grad(XCOORD,iGrad,icg) = grad(XCOORD,iGrad,icg) + wtx*dVar
          grad(YCOORD,iGrad,icg) = grad(YCOORD,iGrad,icg) + wty*dVar
          grad(ZCOORD,iGrad,icg) = grad(ZCOORD,iGrad,icg) + wtz*dVar

          iGrad = iGrad + 1          
        END DO ! iVar
      END DO ! isl
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!                       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!                       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!                       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCellsFast_3D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCellsFast_3D







! ******************************************************************************
!
! Purpose: Compute constrained gradients of any vector or scalar at cell 
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
!   grad        Gradients of variables at cell centers
!
! Notes:
!   1. If the weighting is changed from inverse-distance to none, then the 
!      routine RFLU_ComputeStencilMomentsX in RFLU_ModWeights must also be 
!      adapted. 
!   2. Restricted to linear reconstruction for now.
!   3. Restricted to Dirichlet constraints for now.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCellsConstr(pRegion,iBegVar,iEndVar,iBegGrad, &
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

    INTEGER :: errorFlag,icg,icg2,icl,iCol,ifg,ifl,iGrad,iPatch,isl,iRow,iVar, & 
               nCols,nConstr,nRows,sCount
    INTEGER, DIMENSION(:), ALLOCATABLE :: constrType               
    REAL(RFREAL) :: cwt,dVar,dx,dy,dz,term,varc,gx,gy,gz
    REAL(RFREAL) :: cofg(XCOORD:ZCOORD)
    REAL(RFREAL) :: colMax(4)
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCellsConstr',&
  'RFLU_ModDifferentiationCells.F90')

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::ComputeGradCellsConstr")
#endif

    IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
      CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
    END IF ! iEndVar

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    cwt = pRegion%mixtInput%cReconstCellsWeight

! ******************************************************************************    
!   Loop over cells and compute constrained gradients 
! ******************************************************************************    

    DO icl = 1,pGrid%nCellsConstr
      icg = pGrid%icgConstr(icl)                   

! ==============================================================================
!     Initialize gradients, set cofg, and allocate memory for constraint type
! ==============================================================================

      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = 0.0_RFREAL
        grad(YCOORD,iGrad,icg) = 0.0_RFREAL
        grad(ZCOORD,iGrad,icg) = 0.0_RFREAL
      END DO ! iGrad
          
      cofg(XCOORD) = pGrid%cofg(XCOORD,icg)
      cofg(YCOORD) = pGrid%cofg(YCOORD,icg)
      cofg(ZCOORD) = pGrid%cofg(ZCOORD,icg)            
     
      ALLOCATE(constrType(pGrid%c2cs(icg)%nBFaceMembs),STAT=errorFlag)
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
              
        DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs
          iPatch = pGrid%c2cs(icg)%bFaceMembs(1,isl)        
          ifl    = pGrid%c2cs(icg)%bFaceMembs(2,isl)        
         
          pPatch => pRegion%patches(iPatch)
         
          constrType(isl) = RFLU_GetConstrType(pRegion,pPatch,varInfo(iVar),ifl)
          
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
        
        IF ( nConstr > 0 ) THEN
  
! ------- Allocate temporary memory --------------------------------------------

          nRows = pGrid%c2cs(icg)%nCellMembs + nConstr
          nCols = pRegion%mixtInput%dimens

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
              DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
                dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

                a(isl,1) = term*dx
                a(isl,2) = term*dy
              END DO ! isl  
                            
              iRow = pGrid%c2cs(icg)%nCellMembs              
              
              DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                  iPatch = pGrid%c2cs(icg)%bFaceMembs(1,isl)        
                  ifl    = pGrid%c2cs(icg)%bFaceMembs(2,isl)        
                  
                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl) - cofg(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl) - cofg(YCOORD)

                  term = cwt/SQRT(dx*dx + dy*dy)            

                  iRow = iRow + 1

                  a(iRow,1) = term*dx
                  a(iRow,2) = term*dy
                END IF ! constrType
              END DO ! isl    
              
! --------- Three dimensions              
                                                    
            CASE ( 3 )
              DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
                dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)
                dz = pGrid%cofg(ZCOORD,icg2) - cofg(ZCOORD)                

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                a(isl,1) = term*dx
                a(isl,2) = term*dy
                a(isl,3) = term*dz
              END DO ! isl  
              
              iRow = pGrid%c2cs(icg)%nCellMembs 
              
              DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs              
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN 
                  iPatch = pGrid%c2cs(icg)%bFaceMembs(1,isl)        
                  ifl    = pGrid%c2cs(icg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl) - cofg(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl) - cofg(YCOORD)
                  dz = pPatch%fc(ZCOORD,ifl) - cofg(ZCOORD)                  

                  term = cwt/SQRT(dx*dx + dy*dy + dz*dz)            

                  iRow = iRow + 1

                  a(iRow,1) = term*dx
                  a(iRow,2) = term*dy
                  a(iRow,3) = term*dz
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
            WRITE(*,*) 'ERROR - Singular matrix in RFLU_ComputeGradCellsConstr!'
            STOP
          END IF ! sCount     
! END TEMPORARY
   
! ------- Compute constrained gradients ----------------------------------------

          SELECT CASE ( pRegion%mixtInput%dimens )
          
! --------- Two dimensions          
           
            CASE ( 2 )
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL
        
              DO isl = 1,pGrid%c2cs(icg)%nCellMembs
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
                dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)
 
                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy) 
                           
                dVar = var(iVar,icg2) - var(iVar,icg)

                gx = gx + term*aInv(1,isl)*dVar
                gy = gy + term*aInv(2,isl)*dVar
              END DO ! isl                    
              
              iRow = pGrid%c2cs(icg)%nCellMembs
 
              DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN                   
                  iPatch = pGrid%c2cs(icg)%bFaceMembs(1,isl)        
                  ifl    = pGrid%c2cs(icg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl) - cofg(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl) - cofg(YCOORD)

                  term = cwt/SQRT(dx*dx + dy*dy)            

                  varc = RFLU_GetConstrValue(pRegion,pPatch,varInfo(iVar),ifl)                
                  dVar = varc - var(iVar,icg)

                  iRow = iRow + 1

                  gx = gx + term*aInv(1,iRow)*dVar
                  gy = gy + term*aInv(2,iRow)*dVar
                END IF ! constrType     
              END DO ! isl
                            
              grad(XCOORD,iGrad,icg) = gx
              grad(YCOORD,iGrad,icg) = gy
              grad(ZCOORD,iGrad,icg) = 0.0_RFREAL   
              
! --------- Three dimensions              
                                      
            CASE ( 3 ) 
              gx = 0.0_RFREAL
              gy = 0.0_RFREAL
              gz = 0.0_RFREAL    
    
              DO isl = 1,pGrid%c2cs(icg)%nCellMembs
                icg2 = pGrid%c2cs(icg)%cellMembs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - cofg(XCOORD)
                dy = pGrid%cofg(YCOORD,icg2) - cofg(YCOORD)
                dz = pGrid%cofg(ZCOORD,icg2) - cofg(ZCOORD)                
 
                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

                dVar = var(iVar,icg2) - var(iVar,icg)

                gx = gx + term*aInv(1,isl)*dVar
                gy = gy + term*aInv(2,isl)*dVar
                gz = gz + term*aInv(3,isl)*dVar               
              END DO ! isl                    
              
              iRow = pGrid%c2cs(icg)%nCellMembs

              DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs
                IF ( constrType(isl) == CONSTR_TYPE_DIRICHLET ) THEN               
                  iPatch = pGrid%c2cs(icg)%bFaceMembs(1,isl)        
                  ifl    = pGrid%c2cs(icg)%bFaceMembs(2,isl)        

                  pPatch => pRegion%patches(iPatch)

                  dx = pPatch%fc(XCOORD,ifl) - cofg(XCOORD)
                  dy = pPatch%fc(YCOORD,ifl) - cofg(YCOORD)
                  dz = pPatch%fc(ZCOORD,ifl) - cofg(ZCOORD)                

                  term = cwt/SQRT(dx*dx + dy*dy + dz*dz)   

                  varc = RFLU_GetConstrValue(pRegion,pPatch,varInfo(iVar),ifl)
                  dVar = varc - var(iVar,icg)

                  iRow = iRow + 1

                  gx = gx + term*aInv(1,iRow)*dVar
                  gy = gy + term*aInv(2,iRow)*dVar
                  gz = gz + term*aInv(3,iRow)*dVar
                END IF ! constrType                        
              END DO ! isl
              
              grad(XCOORD,iGrad,icg) = gx
              grad(YCOORD,iGrad,icg) = gy
              grad(ZCOORD,iGrad,icg) = gz                                               
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
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::ComputeGradCellsConstr")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCellsConstr









! ******************************************************************************
!
! Purpose: Compute gradients of any vector or scalar at cell centers. 
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
!   varInfo     Variable information
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_ComputeGradCellsWrapper(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,varInfo,var,grad)

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

    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeGradCellsWrapper',&
  'RFLU_ModDifferentiationCells.F90' )

! ******************************************************************************    
!   Call gradient routines
! ******************************************************************************    

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_ComputeGradCells_1D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                      iEndGrad,var,grad)               
      CASE ( 2 ) 
!        CALL RFLU_ComputeGradCells_2D(pRegion,iBegVar,iEndVar,iBegGrad, &
!                                      iEndGrad,var,grad)
        CALL RFLU_ComputeGradCellsFast_2D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,var,grad)				      
                                         
        IF ( pRegion%grid%nCellsConstr > 0 ) THEN                               
          CALL RFLU_ComputeGradCellsConstr(pRegion,iBegVar,iEndVar,iBegGrad, &
                                           iEndGrad,varInfo,var,grad)
        END IF ! pRegion%grid%nCellsConstr
      CASE ( 3 ) 
!        CALL RFLU_ComputeGradCells_3D(pRegion,iBegVar,iEndVar,iBegGrad, &
!                                      iEndGrad,var,grad)
        CALL RFLU_ComputeGradCellsFast_3D(pRegion,iBegVar,iEndVar,iBegGrad, &
                                          iEndGrad,var,grad)                                      
                                         
        IF ( pRegion%grid%nCellsConstr > 0 ) THEN                               
          CALL RFLU_ComputeGradCellsConstr(pRegion,iBegVar,iEndVar,iBegGrad, &
                                           iEndGrad,varInfo,var,grad)
        END IF ! pRegion%grid%nCellsConstr          
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeGradCellsWrapper





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModDifferentiationCells


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDifferentiationCells.F90,v $
! Revision 1.12  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2007/02/27 13:03:28  haselbac
! Enabled 1d computations
!
! Revision 1.9  2006/04/19 19:42:31  haselbac
! Added tuned gradient routines
!
! Revision 1.8  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.7  2006/04/07 14:48:12  haselbac
! Moved WENO routines into separate module, bug fix for 1D grads
!
! Revision 1.6  2006/03/09 15:05:17  haselbac
! Bug fix: Missing ampersand
!
! Revision 1.5  2006/03/09 14:07:10  haselbac
! Now use generic 1d weight routine
!
! Revision 1.4  2006/01/06 22:10:00  haselbac
! Added wrappers and 1D routines
!
! Revision 1.3  2005/12/25 15:31:33  haselbac
! Added user-specified constraint weight
!
! Revision 1.2  2005/10/27 19:03:25  haselbac
! Separate constrained-gradient routine for performance reasons
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************













