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
! Purpose: Suite of routines related to limiter functions.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModLimiters.F90,v 1.5 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModLimiters

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ComputeLimiterBarthJesp, &
            RFLU_ComputeLimiterVenkat, &
            RFLU_CreateLimiter, &
            RFLU_DestroyLimiter, & 
            RFLU_LimitGradCells, &
            RFLU_LimitGradCellsSimple

        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModLimiters.F90,v $ $Revision: 1.5 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! ******************************************************************************
!
! Purpose: Compute limiters for gradients at cell centers using Barth-
!   Jespersen limiter function.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar	Beginning index in var array
!   iEndVar	Beginning index in var array
!   iBegGrad	Beginning index in grad array
!   iEndGrad	Beginning index in grad array
!   var		Variables to be limited
!   grad	Gradients of variables
!
! Output:
!   lim        	Limiter
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeLimiterBarthJesp(pRegion,iBegVar,iEndVar,iBegGrad, &
                                        iEndGrad,var,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: iBegGrad,iBegVar,iEndGrad,iEndVar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim,var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,errorFlag,icg,ifg,iGrad,iVar
  REAL(RFREAL) :: dx1,dx2,dy1,dy2,dz1,dz2,d1max1,d1max2,d1min1,d1min2,d21, &
                  d22,term,var1,var2,xc,yc,zc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: varMax,varMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLimiterBarthJesp',&
  'RFLU_ModLimiters.F90' )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! iEndVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************    
! Allocate temporary memory
! ******************************************************************************    

  ALLOCATE(varMin(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  ALLOCATE(varMax(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************    
! Initialize
! ******************************************************************************    

  DO icg = 1,pGrid%nCellsTot                       
    DO iVar = iBegVar,iEndVar
      varMax(iVar,icg) = var(iVar,icg) 
      varMin(iVar,icg) = var(iVar,icg)
    END DO ! iVar
    
    DO iGrad = iBegGrad,iEndGrad
      lim(iGrad,icg) = 1.0_RFREAL
    END DO ! iGrad
  END DO ! icg 

! ******************************************************************************    
! Compute varMax and varMin 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)                       

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      varMin(iVar,c1) = MIN(varMin(iVar,c1),var2)
      varMin(iVar,c2) = MIN(varMin(iVar,c2),var1)
      varMax(iVar,c1) = MAX(varMax(iVar,c1),var2)
      varMax(iVar,c2) = MAX(varMax(iVar,c2),var1)
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Loop over faces and compute limiter 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)
  
    dx1 = pGrid%cofg(XCOORD,c1) - xc
    dy1 = pGrid%cofg(YCOORD,c1) - yc
    dz1 = pGrid%cofg(ZCOORD,c1) - zc

    dx2 = pGrid%cofg(XCOORD,c2) - xc
    dy2 = pGrid%cofg(YCOORD,c2) - yc
    dz2 = pGrid%cofg(ZCOORD,c2) - zc

    iGrad = iBegGrad

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      d1min1 = varMin(iVar,c1) - var1
      d1max1 = varMax(iVar,c1) - var1
      d1min2 = varMin(iVar,c2) - var2
      d1max2 = varMax(iVar,c2) - var2

      d21 = grad(XCOORD,iGrad,c1)*dx1 &
          + grad(YCOORD,iGrad,c1)*dy1 &
          + grad(ZCOORD,iGrad,c1)*dz1
      d22 = grad(XCOORD,iGrad,c2)*dx2 &
          + grad(YCOORD,iGrad,c2)*dy2 &
          + grad(ZCOORD,iGrad,c2)*dz2
            
      IF ( d21 > 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1max1/d21)
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ELSEIF ( d21 < 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1min1/d21)
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      END IF ! d21       
                     
      IF ( d22 > 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1max2/d22)
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ELSEIF ( d22 < 0.0_RFREAL ) THEN
        term = MIN(1.0_RFREAL,d1min2/d22)
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      END IF ! d22       

      iGrad = iGrad + 1
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Deallocate temporary memory
! ******************************************************************************    

  DEALLOCATE(varMin,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  DEALLOCATE(varMax,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeLimiterBarthJesp







! ******************************************************************************
!
! Purpose: Compute limiters for gradients at cell centers using Venkatakrishnan
!   limiter function.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar	Beginning index in var array
!   iEndVar	Beginning index in var array
!   iBegGrad	Beginning index in grad array
!   iEndGrad	Beginning index in grad array
!   var		Variables to be limited
!   grad	Gradients of variables
!
! Output:
!   lim        	Limiter
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_ComputeLimiterVenkat(pRegion,iBegVar,iEndVar,iBegGrad, &
                                     iEndGrad,var,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: iBegGrad,iBegVar,iEndGrad,iEndVar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim,var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,errorFlag,icg,ifg,iGrad,iVar
  REAL(RFREAL), PARAMETER :: THRD = 1.0_RFREAL/3.0_RFREAL
  REAL(RFREAL), PARAMETER :: TINY = 1.0E-12_RFREAL 
  REAL(RFREAL) :: denom,ds1,ds2,dx1,dx2,dy1,dy2,dz1,dz2,d1max1,d1max2,d1min1, &
                  d1min2,d21,d22,epsq1,epsq2,numer,term,var1,var2,venkatLimK, &
                  xc,yc,zc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: varMax,varMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLimiterVenkat',&
  'RFLU_ModLimiters.F90' )

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! TEMPORARY
  venkatLimK = 1.0_RFREAL
! END TEMPORARY

! ******************************************************************************    
! Allocate temporary memory
! ******************************************************************************    

  ALLOCATE(varMin(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  ALLOCATE(varMax(iBegVar:iEndVar,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************    
! Initialize
! ******************************************************************************    

  DO icg = 1,pGrid%nCellsTot                       
    DO iVar = iBegVar,iEndVar
      varMax(iVar,icg) = var(iVar,icg) 
      varMin(iVar,icg) = var(iVar,icg)
    END DO ! iVar
    
    DO iGrad = iBegGrad,iEndGrad
      lim(iGrad,icg) = 1.0_RFREAL
    END DO ! iGrad
  END DO ! icg 

! ******************************************************************************    
! Compute varMax and varMin 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)                       

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      varMin(iVar,c1) = MIN(varMin(iVar,c1),var2)
      varMin(iVar,c2) = MIN(varMin(iVar,c2),var1)
      varMax(iVar,c1) = MAX(varMax(iVar,c1),var2)
      varMax(iVar,c2) = MAX(varMax(iVar,c2),var1)
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Loop over faces and compute limiter 
! ******************************************************************************    

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)                       
    c2 = pGrid%f2c(2,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)
  
    dx1 = pGrid%cofg(XCOORD,c1) - xc
    dy1 = pGrid%cofg(YCOORD,c1) - yc
    dz1 = pGrid%cofg(ZCOORD,c1) - zc

    dx2 = pGrid%cofg(XCOORD,c2) - xc
    dy2 = pGrid%cofg(YCOORD,c2) - yc
    dz2 = pGrid%cofg(ZCOORD,c2) - zc

    ds1   = pGrid%vol(c1)**THRD 
    ds2   = pGrid%vol(c2)**THRD 
    epsq1 = (venkatLimK*ds1)*(venkatLimK*ds1)*(venkatLimK*ds1)
    epsq2 = (venkatLimK*ds2)*(venkatLimK*ds2)*(venkatLimK*ds2)

    iGrad = iBegGrad

    DO iVar = iBegVar,iEndVar
      var1 = var(iVar,c1)
      var2 = var(iVar,c2)

      d1min1 = varMin(iVar,c1) - var1
      d1max1 = varMax(iVar,c1) - var1
      d1min2 = varMin(iVar,c2) - var2
      d1max2 = varMax(iVar,c2) - var2

      d21 = grad(XCOORD,iGrad,c1)*dx1 &
          + grad(YCOORD,iGrad,c1)*dy1 &
          + grad(ZCOORD,iGrad,c1)*dz1
      d22 = grad(XCOORD,iGrad,c2)*dx2 &
          + grad(YCOORD,iGrad,c2)*dy2 &
          + grad(ZCOORD,iGrad,c2)*dz2
            
      IF ( d21 > 0.0_RFREAL ) THEN
        numer = (d1max1*d1max1+epsq1)*d21 + 2*d21*d21*d1max1
        denom = d21*(d1max1*d1max1 + 2*d21*d21 + d1max1*d21 + epsq1)
        term  = numer/(denom+TINY)         
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ELSEIF ( d21 < 0.0_RFREAL ) THEN
        numer = (d1min1*d1min1+epsq1)*d21 + 2*d21*d21*d1min1
        denom = d21*(d1min1*d1min1 + 2*d21*d21 + d1min1*d21 + epsq1)
        term  = numer/(denom+TINY) 
        lim(iGrad,c1) = MIN(term,lim(iGrad,c1))
      ENDIF ! d21       
                     
      IF ( d22 > 0.0_RFREAL ) THEN
        numer = (d1max2*d1max2+epsq2)*d22 + 2*d22*d22*d1max2
        denom = d22*(d1max2*d1max2 + 2*d22*d22 + d1max2*d22 + epsq2)
        term  = numer/(denom+TINY) 
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ELSEIF ( d22 < 0.0_RFREAL ) THEN
        numer = (d1min2*d1min2+epsq2)*d22 + 2*d22*d22*d1min2
        denom = d22*(d1min2*d1min2 + 2*d22*d22 + d1min2*d22 + epsq2)
        term  = numer/(denom+TINY) 
        lim(iGrad,c2) = MIN(term,lim(iGrad,c2))
      ENDIF ! d22       

      iGrad = iGrad + 1
    END DO ! iVar
  END DO ! ifg 

! ******************************************************************************    
! Deallocate temporary memory
! ******************************************************************************    

  DEALLOCATE(varMin,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMin')
  END IF ! global%error

  DEALLOCATE(varMax,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varMax')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeLimiterVenkat








! ******************************************************************************
!
! Purpose: Create limiter.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iBegGrad	Beginning index of gradient array
!   iEndGrad	Ending index of gradient array
!
! Output: 
!   lim		Limiter
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_CreateLimiter(pRegion,iBegGrad,iEndGrad,lim)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  CALL RegisterFunction(global,'RFLU_CreateLimiter', &
                        'RFLU_ModLimiters.F90')

! ******************************************************************************
! Create limiter
! ******************************************************************************

  ALLOCATE(lim(iBegGrad:iEndGrad,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'lim')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreateLimiter







! ******************************************************************************
!
! Purpose: Destroy limiter.
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

SUBROUTINE RFLU_DestroyLimiter(pRegion,lim)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyLimiter', &
                        'RFLU_ModLimiters.F90')

! ******************************************************************************
! Destroy limiter
! ******************************************************************************

  DEALLOCATE(lim,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'lim')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyLimiter







! ******************************************************************************
!
! Purpose: Apply limiter to gradient.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iBegGrad	Beginning index of grad array
!   iEndGrad	Ending index of grad array
!   grad	Gradient array
!   lim		Limiter
!
! Output: 
!   grad	Limited gradient array
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_LimitGradCells(pRegion,iBegGrad,iEndGrad,grad,lim) 

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
  REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,iGrad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  CALL RegisterFunction(global,'RFLU_LimitGradCells', &
                        'RFLU_ModLimiters.F90')

! ******************************************************************************
! Apply limiter
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot                       
    DO iGrad = iBegGrad,iEndGrad
      grad(XCOORD,iGrad,icg) = lim(iGrad,icg)*grad(XCOORD,iGrad,icg) 
      grad(YCOORD,iGrad,icg) = lim(iGrad,icg)*grad(YCOORD,iGrad,icg)
      grad(ZCOORD,iGrad,icg) = lim(iGrad,icg)*grad(ZCOORD,iGrad,icg) 
    END DO ! iVar
  END DO ! icg 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_LimitGradCells








! ******************************************************************************
!
! Purpose: Limit gradients of any vector or scalar at cell centers 
!   to ensure positive face quantities.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegVar     Beginning index of data in var
!   iEndVar     Ending index of data in var
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   var         Variables 
!   varInfo     Information on variables
!   grad        Gradients of variables at cell centers
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
SUBROUTINE RFLU_LimitGradCellsSimple(pRegion,iBegVar,iEndVar,iBegGrad, &
                                     iEndGrad,var,varInfo,grad)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iBegVar,iEndVar,iBegGrad,iEndGrad
  INTEGER, DIMENSION(:), POINTER :: varInfo    
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,icg,icl,ict,ifg,ifl,iGrad,iPatch,iVar,nFacesPerCell
  INTEGER, DIMENSION(:,:,:), POINTER :: pC2f
  REAL(RFREAL) :: dx,dy,dz,varc,varf,xc,yc,zc,xfc,yfc,zfc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_LimitGradCellsSimple',&
  'RFLU_ModLimiters.F90' )

  IF ( (iEndVar  - iBegVar) /= (iEndGrad - iBegGrad) ) THEN 
    CALL ErrorStop(global,ERR_GRAD_MISMATCH,__LINE__)
  END IF ! iEndVar

! ******************************************************************************    
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Loop over cells
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    ict = pGrid%cellGlob2Loc(1,icg) ! cell type
    icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

    SELECT CASE ( ict ) 
      CASE ( CELL_TYPE_TET ) 
        nFacesPerCell = 4
        pC2f => pGrid%tet2f              
      CASE ( CELL_TYPE_HEX ) 
        nFacesPerCell = 6
        pC2f => pGrid%hex2f              
      CASE ( CELL_TYPE_PRI ) 
        nFacesPerCell = 5    
        pC2f => pGrid%pri2f                         
      CASE ( CELL_TYPE_PYR ) 
        nFacesPerCell = 5
        pC2f => pGrid%pyr2f                
      CASE DEFAULT  
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict

    xc = pGrid%cofg(XCOORD,icg)
    yc = pGrid%cofg(YCOORD,icg)
    zc = pGrid%cofg(ZCOORD,icg)                    

! ==============================================================================
!   Loop over positive-definite variables
! ==============================================================================     

    iGrad = iBegGrad    

    DO iVar = iBegVar,iEndVar              
      IF ( varInfo(iVar) == VAR_INFO_POS ) THEN 
        varc = var(iVar,icg)

        faceLoop: DO ifl = 1,nFacesPerCell
          iPatch = pC2f(1,ifl,icl)
          ifg    = pC2f(2,ifl,icl)

          IF ( iPatch == 0 ) THEN ! Interior face
            dx = pGrid%fc(XCOORD,ifg) - xc
            dy = pGrid%fc(YCOORD,ifg) - yc
            dz = pGrid%fc(ZCOORD,ifg) - zc                         
          ELSE ! Boundary face
            pPatch => pRegion%patches(iPatch)

            dx = pPatch%fc(XCOORD,ifg) - xc
            dy = pPatch%fc(YCOORD,ifg) - yc
            dz = pPatch%fc(ZCOORD,ifg) - zc                             
          END IF ! pC2f 

          varf = varc + grad(XCOORD,iGrad,icg)*dx &
                      + grad(YCOORD,iGrad,icg)*dy &
                      + grad(ZCOORD,iGrad,icg)*dz

          IF ( varf <= 0.0_RFREAL ) THEN 
            grad(XCOORD,iGrad,icg) = 0.0_RFREAL
            grad(YCOORD,iGrad,icg) = 0.0_RFREAL
            grad(ZCOORD,iGrad,icg) = 0.0_RFREAL

            EXIT faceLoop
          END IF ! varf            

        END DO faceLoop                            
      END IF ! iVar

      iGrad = iGrad + 1        
    END DO ! iVar
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_LimitGradCellsSimple





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModLimiters

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModLimiters.F90,v $
! Revision 1.5  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/27 15:10:55  haselbac
! Many bug fixes and clean-up
!
! Revision 1.2  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.1  2005/07/11 19:33:47  mparmar
! Initial revision
!
! ******************************************************************************












