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
! Purpose: Compute face weights for optimal LES approach.
!
! Description: None.
!
! Input: 
!   region      current region
!
! Output: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: RFLU_ComputeWeightsOLES.F90,v 1.6 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!*******************************************************************************

SUBROUTINE RFLU_ComputeWeightsOLES(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModTools, ONLY: MakeNonZero

  USE RFLU_ModOLES, ONLY: RFLU_GetLPosInvOLES,RFLU_GetQPosInvOLES
  
  USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD  

  IMPLICIT NONE
        
! parameters      
      
  TYPE(t_region) :: region

! local variables

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: a,b,g,hLoc,hOffset,i,ifc,ifcp,ihLoc,ivLoc,j,k,nCells,nCols, & 
             nRows,sCount,vLoc,vOffset       
  REAL(RFREAL) :: lambda,term,wt0,wt1,wt2
  REAL(RFREAL), DIMENSION(:), POINTER :: rhs,vol
  REAL(RFREAL), DIMENSION(:,:), POINTER :: fn,lhs,lhsInv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: int1,int20,int21,int31,int32, & 
                                             int40,int41,int42,int50,int51, & 
                                             int52
  REAL(RFREAL), DIMENSION(:,:,:,:), POINTER :: wtLinOLES
  REAL(RFREAL), DIMENSION(:,:,:,:,:,:), POINTER :: wtQuadOLES
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RFLU_ComputeWeightsOLES',&
  'RFLU_ComputeWeightsOLES.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  int1  => region%grid%int1OLES
  int20 => region%grid%int20OLES
  int21 => region%grid%int21OLES
  int31 => region%grid%int31OLES
  int32 => region%grid%int32OLES
  int40 => region%grid%int40OLES
  int41 => region%grid%int41OLES
  int42 => region%grid%int42OLES  
  int50 => region%grid%int50OLES
  int51 => region%grid%int51OLES
  int52 => region%grid%int52OLES    

  lhs    => region%grid%lhsOLES
  lhsInv => region%grid%lhsInvOLES
  rhs    => region%grid%rhsOLES

  wtLinOLES  => region%grid%wtLinOLES
  wtQuadOLES => region%grid%wtQuadOLES

  fn  => region%grid%fn
  vol => region%grid%vol

! TEMPORARY
!  WRITE(*,*) '### WARNING correcting sign for tangential linear terms ###'
!  WRITE(*,*) '### WARNING skipping weights computation ###'
!  GOTO 1
! END TEMPORARY

! ******************************************************************************
! Compute lambda, used as scaling factor
! ******************************************************************************

! TEMPORARY
  WRITE(*,*) '### WARNING zeroing I3 ###'
! END TEMPORARY

  term = MakeNonZero(global%enerOLES)
  lambda = region%grid%deltaOLES*global%dissOLES/ &  
           (2.0_RFREAL/3.0_RFREAL*term)**(3.0_RFREAL/2.0_RFREAL)   
  
! ******************************************************************************
! Loop over prototype faces
! ******************************************************************************

  DO ifcp = 1,3 ! Loop over prototype faces
    nCells = SIZE(region%grid%fsOLES,1) 
    
! ==============================================================================
!   Compute weights
! ==============================================================================    

    ifc = region%grid%fp2fOLES(ifcp)
               
    wt0 = region%grid%rhoOLES(ifc)/region%grid%deltaOLES
    wt1 = (wt0*lambda)**(2.0_RFREAL/3.0_RFREAL)
    wt2 = (wt0*lambda)**(4.0_RFREAL/3.0_RFREAL)  
        
! DEBUG
    WRITE(*,*) lambda,region%grid%rhoOLES(ifc),region%grid%deltaOLES
    WRITE(*,*) 'weights: ',wt0,wt1,wt2
! END DEBUG        
        
! ==============================================================================
!   Initialize LHS - defensive programming
! ==============================================================================    
    
    lhs(:,:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    
! ==============================================================================    
!   Build LHS
! ==============================================================================    
    
! ------------------------------------------------------------------------------
!   I2
! ------------------------------------------------------------------------------    
                                
    DO ivLoc = 1,3*nCells
      DO ihLoc = 1,3*nCells
        lhs(ivLoc,ihLoc) =     int20(ifcp,ivLoc,ihLoc) & 
                         + wt1*int21(ifcp,ivLoc,ihLoc)                                                   
      END DO ! ihLoc
    END DO ! ivLoc
        
! ------------------------------------------------------------------------------
!   I3 (top right-hand corner)
! ------------------------------------------------------------------------------    
        
    hOffset = 3*nCells
    
    DO ivLoc = 1,3*nCells
      DO ihLoc = 1,9*nCells*nCells
        hLoc = ihLoc + hOffset
             
        lhs(ivLoc,hLoc) = 0.0_RFREAL              
!        lhs(ivLoc,hLoc) = int31(ifcp,ivLoc,ihLoc)        
      END DO ! ihLoc
    END DO ! ivLoc
    
! ------------------------------------------------------------------------------
!   I3 (bottom left-hand corner) - NOTE additional factor involving lambda
! ------------------------------------------------------------------------------    
        
    vOffset = 3*nCells
    
    DO ivLoc = 1,9*nCells*nCells
      DO ihLoc = 1,3*nCells
        vLoc = ivLoc + vOffset
         
        lhs(vLoc,ihLoc) = 0.0_RFREAL            
!        lhs(vLoc,ihLoc) = lambda*lambda*int32(ifcp,ivLoc,ihLoc)
      END DO ! ihLoc
    END DO ! ivLoc    
        
! ------------------------------------------------------------------------------
!   I5
! ------------------------------------------------------------------------------
               
    vOffset = 3*nCells
    hOffset = 3*nCells
        
    DO ivLoc = 1,9*nCells*nCells
      DO ihLoc = 1,9*nCells*nCells
        vLoc = ivLoc + vOffset
        hLoc = ihLoc + hOffset
 
        lhs(vLoc,hLoc) =     int50(ifcp,ivLoc,ihLoc) & 
                       + wt1*int51(ifcp,ivLoc,ihLoc) & 
                       + wt2*int52(ifcp,ivLoc,ihLoc)                                              
      END DO ! ihLoc
    END DO ! ivLoc        
    
! ==============================================================================    
!   Invert LHS using SVD
! ==============================================================================    
     
    nRows = 3*nCells*(1 + 3*nCells)
    nCols = 3*nCells*(1 + 3*nCells)
          
! DEBUG
  WRITE(*,*) 'in RFLU_ComputeWeightsOLES, before SVD'
  WRITE(*,*) MINVAL(int1),MAXVAL(int1)
  WRITE(*,*) MINVAL(int20),MAXVAL(int20)
  WRITE(*,*) MINVAL(int21),MAXVAL(int21)
  WRITE(*,*) MINVAL(int31),MAXVAL(int31)
  WRITE(*,*) MINVAL(int32),MAXVAL(int32) 
  WRITE(*,*) MINVAL(int40),MAXVAL(int40)
  WRITE(*,*) MINVAL(int41),MAXVAL(int41) 
  WRITE(*,*) MINVAL(int42),MAXVAL(int42) 
  WRITE(*,*) MINVAL(int50),MAXVAL(int50)
  WRITE(*,*) MINVAL(int51),MAXVAL(int51) 
  WRITE(*,*) MINVAL(int52),MAXVAL(int52)            
  WRITE(*,*) MINVAL(lhs),MAXVAL(lhs)
! END DEBUG          
          
! DEBUG
!    WRITE(*,*) 'lhs'
!    DO ivLoc = 13,36
!      WRITE(*,'(12(1X,F10.4))') (lhs(ivLoc,ihLoc),ihLoc=1,24)
!    END DO ! ivLoc
! END DEBUG          
          
    CALL RFLU_InvertMatrixSVD(global,nRows,nCols,lhs,lhsInv,sCount)    

! DEBUG
!    IF ( sCount > 0 ) THEN 
      WRITE(*,*) '>>> Matrix singular? <<< ',sCount
!    END IF ! sCount
! END DEBUG

! ==============================================================================    
!   Determine weights for different RHSs
! ==============================================================================    
    
    DO i = 1,3 ! Loop over velocity components
    
! DEBUG
!      WRITE(*,*) '  Component: ',i
! END DEBUG    
    
! ------------------------------------------------------------------------------
!     Initialize RHS
! ------------------------------------------------------------------------------    
    
      rhs(:) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    
! ------------------------------------------------------------------------------
!     Build RHS
! ------------------------------------------------------------------------------    
    
! --- I1 -----------------------------------------------------------------------

      DO ivLoc = 1,3*nCells
        rhs(ivLoc) = int1(i,ifcp,ivLoc)
      END DO ! ivLoc

! --- I4 -----------------------------------------------------------------------

      vOffset = 3*nCells

      DO ivLoc = 1,9*nCells*nCells
        vLoc = ivLoc + vOffset

        rhs(vLoc) =     int40(i,ifcp,ivLoc) & 
                  + wt1*int41(i,ifcp,ivLoc) & 
                  + wt2*int42(i,ifcp,ivLoc)
      END DO ! ivLoc

! DEBUG
!    WRITE(*,*) 'rhs'
!    DO ivLoc = 13,36
!      WRITE(*,'(1(1X,F10.4))') rhs(ivLoc)
!    END DO ! ivLoc
! END DEBUG

! ------------------------------------------------------------------------------
!     Determine, denormalize, and store weights. NOTE that the coefficients 
!     L and Q are defined to act on the product of cell-average quantity and 
!     volume, and because it is easier to compare coefficients which do not 
!     contain the influence of a volume, multiply the coefficients here. If 
!     this happens to get changed at some point in the future, the flux routine
!     also has to change... 
! ------------------------------------------------------------------------------    
    
      rhs = MATMUL(lhsInv,rhs)
    
! DEBUG
!    WRITE(*,*) 'sol'
!    DO ivLoc = 13,36
!      WRITE(*,'(1X,I4,1X,F10.4)') ivLoc,rhs(ivLoc)
!    END DO ! ivLoc
! END DEBUG    
    
! --- Extract linear weights ---------------------------------------------------    
    
      term = MakeNonZero(global%enerOLES)
      term = 3.0_RFREAL*global%dissOLES/(2.0_RFREAL*term)
    
      DO ivLoc = 1,3*nCells
        CALL RFLU_GetLPosInvOLES(ivLoc,j,a)
        wtLinOLES(i,j,a,ifcp) = term*rhs(ivLoc)*vol(a)/fn(XYZMAG,ifc) 
        
        IF ( ABS(wtLinOLES(i,j,a,ifcp)) < 1.0E-12_RFREAL ) THEN 
          wtLinOLES(i,j,a,ifcp) = 0.0_RFREAL
        END IF ! wtLinOLES      
      END DO ! ivLoc
      
! DEBUG
      WRITE(*,*) 'l ',MINVAL(wtLinOLES(i,:,:,ifcp)), & 
                      MAXVAL(wtLinOLES(i,:,:,ifcp))
! END DEBUG       
      
! --- Extract quadratic weights ------------------------------------------------      
            
      term = 1.0_RFREAL/region%grid%deltaOLES**4      
                     
      vOffset = 3*nCells      
                     
      DO ivLoc = 1,9*nCells*nCells
        vLoc = ivLoc + vOffset
        
        CALL RFLU_GetQPosInvOLES(ivLoc,nCells,j,k,b,g) 
        wtQuadOLES(i,j,k,b,g,ifcp) = term*vol(b)*vol(g)*rhs(vLoc)/fn(XYZMAG,ifc)
           
        IF ( ABS(wtQuadOLES(i,j,k,b,g,ifcp)) < 1.0E-12_RFREAL ) THEN 
          wtQuadOLES(i,j,k,b,g,ifcp) = 0.0_RFREAL
        END IF ! wtQuadOLES     
      END DO ! ivLoc
            
! DEBUG
      WRITE(*,*) 'q ',MINVAL(wtQuadOLES(i,:,:,:,:,ifcp)), & 
                      MAXVAL(wtQuadOLES(i,:,:,:,:,ifcp))
! END DEBUG              
            
! DEBUG
!      WRITE(*,*) 'Local face:  ',ifcp
!      WRITE(*,*) '  Global face: ',ifc 
!      WRITE(*,*) '  Component:   ',i
!      WRITE(*,'(1X,A,3(1X,E13.6))') '  Normal:      ',fn(1:3,ifc)
!      WRITE(*,*) '  Stencil:     ',region%grid%fsOLES(1:nCells,ifc)
!
!      WRITE(*,*) '  Linear part: '
!      DO ivLoc = 1,3*nCells
!        CALL RFLU_GetLPosInvOLES(ivLoc,j,a)
!        
!        WRITE(*,*) '    Counter:     ',ivLoc
!        WRITE(*,*) '    Indices:     ',j,a
!!        WRITE(*,*) '    Global cell: ',region%grid%fsOLES(a,ifc)
!        WRITE(*,*) '    Weight:      ',wtLinOLES(i,j,a,ifcp)                
!      END DO ! ivLoc
!      
!      WRITE(*,*) '  Quadratic part: '
!      DO ivLoc = 1,9*nCells*nCells
!        CALL RFLU_GetQPosInvOLES(ivLoc,nCells,j,k,b,g)
!        
!        WRITE(*,*) '    Counter:      ',ivLoc
!        WRITE(*,*) '    Indices:      ',j,k,b,g
!!        WRITE(*,*) '    Global cells: ',region%grid%fsOLES(b,ifc),region%grid%fsOLES(g,ifc)
!        WRITE(*,*) '    Weights:      ',wtQuadOLES(i,j,k,b,g,ifcp)                
!      END DO ! ivLoc        
! END DEBUG       
                        
    END DO ! i
  END DO ! ifcp

! TEMPORARY
1 CONTINUE

! Jacobs weights
!  GOTO 2
  GOTO 3
  WRITE(*,*) '### WARNING using Jacobs weights ###'

  wtLinOLES(:,:,:,:) = 0.0_RFREAL
  
  wtLinOLES(1,1,1,1) =  0.357_RFREAL
  wtLinOLES(1,1,2,1) = -0.357_RFREAL    
  wtLinOLES(2,2,1,1) =  0.156_RFREAL
  wtLinOLES(2,2,2,1) = -0.156_RFREAL
  wtLinOLES(3,3,1,1) =  0.156_RFREAL
  wtLinOLES(3,3,2,1) = -0.156_RFREAL

  wtLinOLES(1,1,1,2) =  0.156_RFREAL
  wtLinOLES(1,1,2,2) = -0.156_RFREAL
  wtLinOLES(2,2,1,2) =  0.357_RFREAL
  wtLinOLES(2,2,2,2) = -0.357_RFREAL 
  wtLinOLES(3,3,1,2) =  0.156_RFREAL
  wtLinOLES(3,3,2,2) = -0.156_RFREAL

  wtLinOLES(1,1,1,3) =  0.156_RFREAL
  wtLinOLES(1,1,2,3) = -0.156_RFREAL
  wtLinOLES(2,2,1,3) =  0.156_RFREAL
  wtLinOLES(2,2,2,3) = -0.156_RFREAL 
  wtLinOLES(3,3,1,3) =  0.357_RFREAL
  wtLinOLES(3,3,2,3) = -0.357_RFREAL

  wtQuadOLES(:,:,:,:,:,:) = 0.0_RFREAL

  wtQuadOLES(1,1,1,1,1,1) =  0.388_RFREAL
  wtQuadOLES(1,1,1,1,2,1) =  0.135_RFREAL  
  wtQuadOLES(1,1,1,2,1,1) =  0.135_RFREAL 
  wtQuadOLES(1,1,1,2,2,1) =  0.388_RFREAL
  
  wtQuadOLES(2,1,2,1,1,1) =  0.180_RFREAL 
  wtQuadOLES(2,2,1,1,1,1) =  0.180_RFREAL
  wtQuadOLES(2,1,2,1,2,1) =  0.087_RFREAL 
  wtQuadOLES(2,2,1,1,2,1) =  0.087_RFREAL   
  wtQuadOLES(2,1,2,2,1,1) =  0.087_RFREAL 
  wtQuadOLES(2,2,1,2,1,1) =  0.087_RFREAL 
  wtQuadOLES(2,1,2,2,2,1) =  0.180_RFREAL 
  wtQuadOLES(2,2,1,2,2,1) =  0.180_RFREAL

  wtQuadOLES(3,1,3,1,1,1) =  0.180_RFREAL 
  wtQuadOLES(3,3,1,1,1,1) =  0.180_RFREAL
  wtQuadOLES(3,1,3,1,2,1) =  0.087_RFREAL 
  wtQuadOLES(3,3,1,1,2,1) =  0.087_RFREAL   
  wtQuadOLES(3,1,3,2,1,1) =  0.087_RFREAL 
  wtQuadOLES(3,3,1,2,1,1) =  0.087_RFREAL 
  wtQuadOLES(3,1,3,2,2,1) =  0.180_RFREAL 
  wtQuadOLES(3,3,1,2,2,1) =  0.180_RFREAL
  
  wtQuadOLES(2,2,2,1,1,2) =  0.388_RFREAL
  wtQuadOLES(2,2,2,1,2,2) =  0.135_RFREAL  
  wtQuadOLES(2,2,2,2,1,2) =  0.135_RFREAL 
  wtQuadOLES(2,2,2,2,2,2) =  0.388_RFREAL  
  
  wtQuadOLES(1,1,2,1,1,2) =  0.180_RFREAL 
  wtQuadOLES(1,2,1,1,1,2) =  0.180_RFREAL
  wtQuadOLES(1,1,2,1,2,2) =  0.087_RFREAL 
  wtQuadOLES(1,2,1,1,2,2) =  0.087_RFREAL   
  wtQuadOLES(1,1,2,2,1,2) =  0.087_RFREAL 
  wtQuadOLES(1,2,1,2,1,2) =  0.087_RFREAL 
  wtQuadOLES(1,1,2,2,2,2) =  0.180_RFREAL 
  wtQuadOLES(1,2,1,2,2,2) =  0.180_RFREAL

  wtQuadOLES(3,2,3,1,1,2) =  0.180_RFREAL 
  wtQuadOLES(3,3,2,1,1,2) =  0.180_RFREAL
  wtQuadOLES(3,2,3,1,2,2) =  0.087_RFREAL 
  wtQuadOLES(3,3,2,1,2,2) =  0.087_RFREAL   
  wtQuadOLES(3,2,3,2,1,2) =  0.087_RFREAL 
  wtQuadOLES(3,3,2,2,1,2) =  0.087_RFREAL 
  wtQuadOLES(3,2,3,2,2,2) =  0.180_RFREAL 
  wtQuadOLES(3,3,2,2,2,2) =  0.180_RFREAL

  wtQuadOLES(3,3,3,1,1,3) =  0.388_RFREAL
  wtQuadOLES(3,3,3,1,2,3) =  0.135_RFREAL  
  wtQuadOLES(3,3,3,2,1,3) =  0.135_RFREAL 
  wtQuadOLES(3,3,3,2,2,3) =  0.388_RFREAL  
  
  wtQuadOLES(1,1,3,1,1,3) =  0.180_RFREAL 
  wtQuadOLES(1,3,1,1,1,3) =  0.180_RFREAL
  wtQuadOLES(1,1,3,1,2,3) =  0.087_RFREAL 
  wtQuadOLES(1,3,1,1,2,3) =  0.087_RFREAL   
  wtQuadOLES(1,1,3,2,1,3) =  0.087_RFREAL 
  wtQuadOLES(1,3,1,2,1,3) =  0.087_RFREAL 
  wtQuadOLES(1,1,3,2,2,3) =  0.180_RFREAL 
  wtQuadOLES(1,3,1,2,2,3) =  0.180_RFREAL

  wtQuadOLES(2,2,3,1,1,3) =  0.180_RFREAL 
  wtQuadOLES(2,3,2,1,1,3) =  0.180_RFREAL
  wtQuadOLES(2,2,3,1,2,3) =  0.087_RFREAL 
  wtQuadOLES(2,3,2,1,2,3) =  0.087_RFREAL   
  wtQuadOLES(2,2,3,2,1,3) =  0.087_RFREAL 
  wtQuadOLES(2,3,2,2,1,3) =  0.087_RFREAL 
  wtQuadOLES(2,2,3,2,2,3) =  0.180_RFREAL 
  wtQuadOLES(2,3,2,2,2,3) =  0.180_RFREAL  
  
! Central flux weights
  GOTO 3
2 CONTINUE
  WRITE(*,*) '### WARNING using central weights ###'
  wtLinOLES(:,:,:,:) = 0.0_RFREAL
 
  wtLinOLES(1,1,1,1) =  0.5_RFREAL
  wtLinOLES(1,1,2,1) = -0.5_RFREAL    
  wtLinOLES(2,2,1,1) =  0.5_RFREAL
  wtLinOLES(2,2,2,1) = -0.5_RFREAL
  wtLinOLES(3,3,1,1) =  0.5_RFREAL
  wtLinOLES(3,3,2,1) = -0.5_RFREAL

  wtLinOLES(1,1,1,2) =  0.5_RFREAL
  wtLinOLES(1,1,2,2) = -0.5_RFREAL
  wtLinOLES(2,2,1,2) =  0.5_RFREAL
  wtLinOLES(2,2,2,2) = -0.5_RFREAL 
  wtLinOLES(3,3,1,2) =  0.5_RFREAL
  wtLinOLES(3,3,2,2) = -0.5_RFREAL

  wtLinOLES(1,1,1,3) =  0.5_RFREAL
  wtLinOLES(1,1,2,3) = -0.5_RFREAL
  wtLinOLES(2,2,1,3) =  0.5_RFREAL
  wtLinOLES(2,2,2,3) = -0.5_RFREAL 
  wtLinOLES(3,3,1,3) =  0.5_RFREAL
  wtLinOLES(3,3,2,3) = -0.5_RFREAL 
 
  wtQuadOLES(:,:,:,:,:,:) = 0.0_RFREAL

! x-face

  wtQuadOLES(1,1,1,1,1,1) = 0.250_RFREAL 
  wtQuadOLES(1,1,1,1,2,1) = 0.250_RFREAL 
  wtQuadOLES(1,1,1,2,2,1) = 0.250_RFREAL 
  wtQuadOLES(1,1,1,2,1,1) = 0.250_RFREAL
   
  wtQuadOLES(2,1,2,1,1,1) = 0.125_RFREAL  
  wtQuadOLES(2,1,2,1,2,1) = 0.125_RFREAL
  wtQuadOLES(2,2,1,1,1,1) = 0.125_RFREAL            
  wtQuadOLES(2,2,1,2,1,1) = 0.125_RFREAL
  wtQuadOLES(2,1,2,2,1,1) = 0.125_RFREAL          
  wtQuadOLES(2,2,1,1,2,1) = 0.125_RFREAL 
  wtQuadOLES(2,1,2,2,2,1) = 0.125_RFREAL
  wtQuadOLES(2,2,1,2,2,1) = 0.125_RFREAL    
  
  wtQuadOLES(3,1,3,1,1,1) = 0.125_RFREAL  
  wtQuadOLES(3,1,3,1,2,1) = 0.125_RFREAL
  wtQuadOLES(3,3,1,1,1,1) = 0.125_RFREAL            
  wtQuadOLES(3,3,1,2,1,1) = 0.125_RFREAL
  wtQuadOLES(3,1,3,2,1,1) = 0.125_RFREAL          
  wtQuadOLES(3,3,1,1,2,1) = 0.125_RFREAL 
  wtQuadOLES(3,1,3,2,2,1) = 0.125_RFREAL
  wtQuadOLES(3,3,1,2,2,1) = 0.125_RFREAL 

! y-face

  wtQuadOLES(2,2,2,1,1,2) = 0.250_RFREAL 
  wtQuadOLES(2,2,2,1,2,2) = 0.250_RFREAL 
  wtQuadOLES(2,2,2,2,2,2) = 0.250_RFREAL 
  wtQuadOLES(2,2,2,2,1,2) = 0.250_RFREAL   
   
  wtQuadOLES(1,1,2,1,1,2) = 0.125_RFREAL  
  wtQuadOLES(1,1,2,1,2,2) = 0.125_RFREAL
  wtQuadOLES(1,2,1,1,1,2) = 0.125_RFREAL            
  wtQuadOLES(1,2,1,2,1,2) = 0.125_RFREAL
  wtQuadOLES(1,1,2,2,1,2) = 0.125_RFREAL          
  wtQuadOLES(1,2,1,1,2,2) = 0.125_RFREAL 
  wtQuadOLES(1,1,2,2,2,2) = 0.125_RFREAL
  wtQuadOLES(1,2,1,2,2,2) = 0.125_RFREAL    
  
  wtQuadOLES(3,2,3,1,1,2) = 0.125_RFREAL  
  wtQuadOLES(3,2,3,1,2,2) = 0.125_RFREAL
  wtQuadOLES(3,3,2,1,1,2) = 0.125_RFREAL            
  wtQuadOLES(3,3,2,2,1,2) = 0.125_RFREAL
  wtQuadOLES(3,2,3,2,1,2) = 0.125_RFREAL          
  wtQuadOLES(3,3,2,1,2,2) = 0.125_RFREAL 
  wtQuadOLES(3,2,3,2,2,2) = 0.125_RFREAL
  wtQuadOLES(3,3,2,2,2,2) = 0.125_RFREAL   
   
! z-face
   
  wtQuadOLES(3,3,3,1,1,3) = 0.250_RFREAL 
  wtQuadOLES(3,3,3,1,2,3) = 0.250_RFREAL 
  wtQuadOLES(3,3,3,2,2,3) = 0.250_RFREAL 
  wtQuadOLES(3,3,3,2,1,3) = 0.250_RFREAL   
   
  wtQuadOLES(1,1,3,1,1,3) = 0.125_RFREAL  
  wtQuadOLES(1,1,3,1,2,3) = 0.125_RFREAL
  wtQuadOLES(1,3,1,1,1,3) = 0.125_RFREAL            
  wtQuadOLES(1,3,1,2,1,3) = 0.125_RFREAL
  wtQuadOLES(1,1,3,2,1,3) = 0.125_RFREAL          
  wtQuadOLES(1,3,1,1,2,3) = 0.125_RFREAL 
  wtQuadOLES(1,1,3,2,2,3) = 0.125_RFREAL
  wtQuadOLES(1,3,1,2,2,3) = 0.125_RFREAL    
  
  wtQuadOLES(2,2,3,1,1,3) = 0.125_RFREAL  
  wtQuadOLES(2,2,3,1,2,3) = 0.125_RFREAL
  wtQuadOLES(2,3,2,1,1,3) = 0.125_RFREAL            
  wtQuadOLES(2,3,2,2,1,3) = 0.125_RFREAL
  wtQuadOLES(2,2,3,2,1,3) = 0.125_RFREAL          
  wtQuadOLES(2,3,2,1,2,3) = 0.125_RFREAL 
  wtQuadOLES(2,2,3,2,2,3) = 0.125_RFREAL
  wtQuadOLES(2,3,2,2,2,3) = 0.125_RFREAL 
3 CONTINUE  
! END TEMPORARY


! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeWeightsOLES

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeWeightsOLES.F90,v $
! Revision 1.6  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2003/01/28 14:33:14  haselbac
! Use parameters in fn
!
! Revision 1.3  2002/09/09 15:45:10  haselbac
! global now under regions, several bug fixes
!
!*******************************************************************************







