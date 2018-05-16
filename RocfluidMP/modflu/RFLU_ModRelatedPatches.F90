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
! Purpose: Suite of routines for transforming data on related patches.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModRelatedPatches.F90,v 1.5 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRelatedPatches

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_RELP_TransformVector, & 
            RFLU_RELP_TransformWrapper
            
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRelatedPatches.F90,v $ $Revision: 1.5 $' 
                      
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Wrapper for transforming data for virtual cells associated with
!   periodic patches.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch      Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RELP_PeriodicWrapper(pRegion,pPatch)

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
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl
    REAL(RFREAL) :: ct,ex,ey,ez,st,theta
    REAL(RFREAL) :: v(3),vr(3)
    REAL(RFREAL) :: rotmat(3,3)
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pVar    
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RELP_PeriodicWrapper',&
  'RFLU_ModRelatedPatches.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid
        
! ******************************************************************************
!   Check boundary condition - defensive coding
! ******************************************************************************

    IF ( pPatch%bcType /= BC_PERIODIC ) THEN 
      CALL ErrorStop(global,ERR_BC_INVALID,__LINE__)
    END IF ! pPatch%bcType

! ******************************************************************************
!   Define rotation matrix
! ******************************************************************************
    
    theta = pPatch%angleRelated

    ct = COS(theta)
    st = SIN(theta)

    SELECT CASE ( pPatch%axisRelated ) 
      CASE ( 1 ) 
        ex = 1.0_RFREAL
        ey = 0.0_RFREAL
        ez = 0.0_RFREAL
      CASE ( 2 ) 
        ex = 0.0_RFREAL
        ey = 1.0_RFREAL
        ez = 0.0_RFREAL
      CASE ( 3 ) 
        ex = 0.0_RFREAL
        ey = 0.0_RFREAL
        ez = 1.0_RFREAL
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%axisRelated
    
    rotmat(1,1) = ct + (1.0_RFREAL-ct)*ex*ex
    rotmat(1,2) =      (1.0_RFREAL-ct)*ex*ey - st*ez
    rotmat(1,3) =      (1.0_RFREAL-ct)*ex*ez + st*ey
    
    rotmat(2,1) =      (1.0_RFREAL-ct)*ey*ex + st*ez
    rotmat(2,2) = ct + (1.0_RFREAL-ct)*ey*ey
    rotmat(2,3) =      (1.0_RFREAL-ct)*ey*ez - st*ex
    
    rotmat(3,1) =      (1.0_RFREAL-ct)*ez*ex - st*ey
    rotmat(3,2) =      (1.0_RFREAL-ct)*ez*ey + st*ex
    rotmat(3,3) = ct + (1.0_RFREAL-ct)*ez*ez        

! ******************************************************************************
!   Mixture
! ******************************************************************************
    
    pVar => pRegion%mixt%cv
    
    DO icl = 1,pPatch%nBCellsVirt
      icg = pPatch%bvc(icl)
      
      v(1) = pVar(CV_MIXT_XMOM,icg)
      v(2) = pVar(CV_MIXT_YMOM,icg)
      v(3) = pVar(CV_MIXT_ZMOM,icg)      
      
      vr = MATMUL(rotmat,v)
            
      pVar(CV_MIXT_XMOM,icg) = vr(1)
      pVar(CV_MIXT_YMOM,icg) = vr(2)
      pVar(CV_MIXT_ZMOM,icg) = vr(3)   
    END DO ! icl

! ******************************************************************************
!   Physical modules
! ******************************************************************************

! TO DO 
!
! END TO DO 

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RELP_PeriodicWrapper






! ******************************************************************************
!
! Purpose: Wrapper for transforming data for virtual cells associated with
!   symmetry patches.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch      Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RELP_SymmetryWrapper(pRegion,pPatch)

    USE ModInterfaces, ONLY: ReflectVector

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
!   Local variables
! ==============================================================================

    INTEGER :: icg,icl
    REAL(RFREAL) :: nx,ny,nz,vx,vy,vz
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pVar
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RELP_SymmetryWrapper',&
  'RFLU_ModRelatedPatches.F90')
        
! ******************************************************************************
!   Check boundary condition - defensive coding
! ******************************************************************************

    IF ( pPatch%bcType /= BC_SYMMETRY ) THEN 
      CALL ErrorStop(global,ERR_BC_INVALID,__LINE__)
    END IF ! pPatch%bcType

! ******************************************************************************
!   Get patch normal
! ******************************************************************************
    
    nx = pPatch%pn(XCOORD)
    ny = pPatch%pn(YCOORD)
    nz = pPatch%pn(ZCOORD)

! ******************************************************************************
!   Mixture
! ******************************************************************************
    
    pVar => pRegion%mixt%cv
    
    DO icl = 1,pPatch%nBCellsVirt
      icg = pPatch%bvc(icl)
      
      vx = pVar(CV_MIXT_XMOM,icg)
      vy = pVar(CV_MIXT_YMOM,icg)
      vz = pVar(CV_MIXT_ZMOM,icg)      
      
      CALL ReflectVector(nx,ny,nz,vx,vy,vz)
      
      pVar(CV_MIXT_XMOM,icg) = vx
      pVar(CV_MIXT_YMOM,icg) = vy
      pVar(CV_MIXT_ZMOM,icg) = vz   
    END DO ! icl

! ******************************************************************************
!   Physical modules
! ******************************************************************************

! TO DO 
!
! END TO DO 

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RELP_SymmetryWrapper







! ******************************************************************************
!
! Purpose: Transform vector.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch      Pointer to patch
!   vx          x-component of vector
!   vy          y-component of vector
!   vz          z-component of vector
!
! Output: 
!   vx          x-component of transformed vector
!   vy          y-component of transformed vector
!   vz          z-component of transformed vector
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_RELP_TransformVector(pRegion,pPatch,vx,vy,vz)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(INOUT) :: vx,vy,vz 
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    REAL(RFREAL) :: v(4),vt(4)    
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RELP_TransformVector',&
  'RFLU_ModRelatedPatches.F90')

! ******************************************************************************
!   Check boundary condition - defensive coding
! ******************************************************************************

    IF ( pPatch%bcType /= BC_PERIODIC ) THEN 
      CALL ErrorStop(global,ERR_BC_INVALID,__LINE__)
    END IF ! pPatch%bcType

! ******************************************************************************
!   Transform vector
! ******************************************************************************
    
    v(1) = vx
    v(2) = vy
    v(3) = vz
    v(4) = 1.0_RFREAL
           
    vt = MATMUL(pPatch%tm,v)
      
    vx = vt(1)
    vy = vt(2)
    vz = vt(3)        

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RELP_TransformVector






! ******************************************************************************
!
! Purpose: Wrapper for transforming data for virtual cells associated with
!   periodic patches.
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

  SUBROUTINE RFLU_RELP_TransformWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
    TYPE(t_patch), POINTER :: pPatch
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_RELP_TransformWrapper',&
  'RFLU_ModRelatedPatches.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
    
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        SELECT CASE ( pPatch%bcType ) 
          CASE ( BC_PERIODIC ) 
            CALL RFLU_RELP_PeriodicWrapper(pRegion,pPatch)
          CASE ( BC_SYMMETRY ) 
            CALL RFLU_RELP_SymmetryWrapper(pRegion,pPatch)          
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType      
      END IF ! pPatch%nBCellsVirt
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_RELP_TransformWrapper









! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModRelatedPatches


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRelatedPatches.F90,v $
! Revision 1.5  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 16:09:51  fnajjar
! Added routine to transform vector
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2006/03/25 21:38:54  haselbac
! Initial revision
!
! ******************************************************************************










