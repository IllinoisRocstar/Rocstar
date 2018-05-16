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
! Purpose: Compute second-order accurate discretization of scalar inviscid flux
!   through boundary faces.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to data of current region
!   pPatch              Pointer to data of current patch
!   nVarScal            Number of scalars
!   cvScal              Vector of conserved scalar variables
!   gradCellScal        Cell gradients of scalar variables
!   valScal             Boundary values of scalar variables     
!
! Output: 
!   resScal             Residual of scalar variables
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ScalarSecondPatch.F90,v 1.8 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ScalarSecondPatch(pRegion,pPatch,nVarScal,cvScal, & 
                                  gradCellScal,valScal,resScal)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_bcvalues,t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
    
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
  REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradCellScal  
  TYPE(t_bcvalues) :: valScal
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,bcType,distScal,ifc,iVarScal
  REAL(RFREAL) :: dx,dy,dz,flx,mf,sl
  REAL(RFREAL), DIMENSION(:), POINTER :: pMfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: rhs 
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarSecondPatch.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarSecondPatch',&
  'RFLU_ScalarSecondPatch.F90')

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************
 
  IF ( pRegion%mixtInput%indMfMixt /= 1 ) THEN 
    CALL ErrorStop(global,ERR_INDMFMIXT_INVALID,__LINE__)
  END IF ! pRegion%mixtInput%indMfMixt
  
! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
      
  pGrid   => pRegion%grid     
  pMfMixt => pPatch%mfMixt
  
  bcType   = pPatch%bcType  
  distScal = valScal%distrib  

! ******************************************************************************
! Select boundary type
! ******************************************************************************

  SELECT CASE ( bcType )

! ==============================================================================  
!   Inflow
! ==============================================================================

    CASE ( BC_INFLOW_TOTANG,BC_INFLOW_VELTEMP ) 
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)
                    
        mf = pMfMixt(ifc)
                                                                        
        DO iVarScal = 1,nVarScal                                      
          flx = mf*valScal%vals(iVarScal,distScal*ifc)
          
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal
      END DO ! ifc
            
! ==============================================================================  
!   Outflow
! ==============================================================================

    CASE ( BC_OUTFLOW )           
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)
        
        dx = pPatch%fc(XCOORD,ifc) - pGrid%cofg(XCOORD,c1)
        dy = pPatch%fc(YCOORD,ifc) - pGrid%cofg(YCOORD,c1)
        dz = pPatch%fc(ZCOORD,ifc) - pGrid%cofg(ZCOORD,c1)            

        mf = pMfMixt(ifc)   
                                                                                                        
        DO iVarScal = 1,nVarScal  
          sl = cvScal(iVarScal,c1)
                      
          sl = sl + gradCellScal(XCOORD,iVarScal,c1)*dx & 
                  + gradCellScal(YCOORD,iVarScal,c1)*dy & 
                  + gradCellScal(ZCOORD,iVarScal,c1)*dz            
                      
          flx = mf*sl
                                  
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal 
      END DO ! ifc

! ==============================================================================  
!   Slip wall 
! ==============================================================================
    
    CASE ( BC_SLIPWALL ) 

! ==============================================================================  
!   No-slip wall
! ==============================================================================
 
    CASE ( BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP ) 

! ==============================================================================  
!   Farfield
! ==============================================================================

    CASE ( BC_FARFIELD )
      DO ifc = 1,pPatch%nBFaces            
        c1 = pPatch%bf2c(ifc)      

        mf = pMfMixt(ifc)   

        IF ( mf > 0.0_RFREAL ) THEN ! Outflow
          dx = pPatch%fc(XCOORD,ifc) - pGrid%cofg(XCOORD,c1)
          dy = pPatch%fc(YCOORD,ifc) - pGrid%cofg(YCOORD,c1)
          dz = pPatch%fc(ZCOORD,ifc) - pGrid%cofg(ZCOORD,c1)            
                  
          DO iVarScal = 1,nVarScal        
            sl = cvScal(iVarScal,c1)

            sl = sl + gradCellScal(XCOORD,iVarScal,c1)*dx & 
                    + gradCellScal(YCOORD,iVarScal,c1)*dy & 
                    + gradCellScal(ZCOORD,iVarScal,c1)*dz           

            flx = mf*sl 
                        
            resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
          END DO ! iVarScal                      
        ELSE ! Inflow
          DO iVarScal = 1,nVarScal 
            flx = mf*valScal%vals(iVarScal,distScal*ifc)
                      
            resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
          END DO ! iVarScal
        END IF ! mf           
      END DO ! ifc      
    
! ==============================================================================  
!   Injection
! ==============================================================================

    CASE ( BC_INJECTION )
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)
             
        mf = pMfMixt(ifc)               
                                        
        DO iVarScal = 1,nVarScal                
          flx = mf*valScal%vals(iVarScal,distScal*ifc)
                                  
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal                           
      END DO ! ifc

! ==============================================================================
!   Boundaries for which fluxes must not or need not be computed
! ==============================================================================

    CASE ( BC_PERIODIC, &
           BC_SYMMETRY, & 
           BC_VIRTUAL )
           
! ==============================================================================  
!   Default
! ==============================================================================

    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarSecondPatch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarSecondPatch.F90,v $
! Revision 1.8  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.5  2006/03/26 20:28:42  haselbac
! Removed unnecessary BC_RANGE
!
! Revision 1.4  2006/03/25 21:43:47  haselbac
! Added CASEs for sype patches
!
! Revision 1.3  2005/11/10 02:03:16  haselbac
! Added virtual boundary, cleaned up CASE statements
!
! Revision 1.2  2005/04/20 14:40:15  haselbac
! Removed CHECK_UNIFLOW code section, cosmetics
!
! Revision 1.1  2004/01/29 22:56:14  haselbac
! Initial revision
!
! ******************************************************************************







