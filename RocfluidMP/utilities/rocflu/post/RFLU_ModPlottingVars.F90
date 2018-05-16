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
! Purpose: Collection of routines for handling plotting variables.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModPlottingVars.F90,v 1.15 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPlottingVars

  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region

  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsWrapper

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

  PRIVATE

! ==============================================================================
! Private data
! ==============================================================================
 
  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPlottingVars.F90,v $ $Revision: 1.15 $'

! ==============================================================================
! Public data
! ==============================================================================


! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_BuildPlottingVarMaps, &
            RFLU_ComputePlottingVarsWrapper, &  
            RFLU_CountPlottingVars, & 
            RFLU_CreatePlottingVarMaps, &
            RFLU_CreatePlottingVars, & 
            RFLU_CreatePlottingVarsVert, & 
            RFLU_DecideComputePlottingVars, &
            RFLU_DestroyPlottingVarMaps, & 
            RFLU_DestroyPlottingVars, & 
            RFLU_DestroyPlottingVarsVert, & 
            RFLU_PrintPlottingVarsInfo

! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Build plotting variable maps.
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

SUBROUTINE RFLU_BuildPlottingVarMaps(pRegion)

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

  INTEGER :: errorFlag,iPv,iPv2,nPv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_BuildPlottingVarMaps', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building plotting variable maps...'
  END IF ! global%verbLevel

! ******************************************************************************
! Build map
! ******************************************************************************

  nPv = 0

  DO iPv = 1,PV_XXXX_NVAR
    iPv2 = pRegion%plot%pv2pvi(iPv)
  
    IF ( iPv2 /= CRAZY_VALUE_INT ) THEN
      nPv = nPv + 1
     
      pRegion%plot%pvi2pv(iPv2) = iPv
    END IF ! iPv2 
  END DO ! iPv

  IF ( nPv /= pRegion%plot%nPv ) THEN 
! TEMPORARY  
    WRITE(*,*) 'ERROR!!!'
    STOP
! END TEMPORARY
  END IF ! nPv

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building plotting variable maps done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BuildPlottingVarMaps 
 
 
 
 

! ******************************************************************************
!
! Purpose: Set plotting variable names.
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

SUBROUTINE RFLU_BuildPlottingVarNames(pRegion)

#ifdef PLAG
  USE PLAG_ModParameters
#endif

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

  INTEGER :: errorFlag,iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_BuildPlottingVarNames', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Set names of plotting variables
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================
        
! ------------------------------------------------------------------------------
! Discontinuity variables
! ------------------------------------------------------------------------------
  
  pRegion%plot%pvNameShort(PV_MIXT_SCHL) = "schl"
  pRegion%plot%pvNameShort(PV_MIXT_SHAD) = "shad"
  pRegion%plot%pvNameShort(PV_MIXT_INTF) = "intf"   
  
  pRegion%plot%pvNameLong(PV_MIXT_SCHL) = "Schlieren variable"  
  pRegion%plot%pvNameLong(PV_MIXT_SHAD) = "Shadowgraph variable"
  pRegion%plot%pvNameLong(PV_MIXT_INTF) = "Interferogram variable"     

! ------------------------------------------------------------------------------
! Vorticity
! ------------------------------------------------------------------------------

  pRegion%plot%pvNameShort(PV_MIXT_XVOR) = "vorx"
  pRegion%plot%pvNameShort(PV_MIXT_YVOR) = "vory"
  pRegion%plot%pvNameShort(PV_MIXT_ZVOR) = "vorz"

  pRegion%plot%pvNameLong(PV_MIXT_XVOR) = "Vorticity x-component"  
  pRegion%plot%pvNameLong(PV_MIXT_YVOR) = "Vorticity y-component"  
  pRegion%plot%pvNameLong(PV_MIXT_ZVOR) = "Vorticity z-component"  

! ------------------------------------------------------------------------------
! Vortex core variables
! ------------------------------------------------------------------------------

  pRegion%plot%pvNameShort(PV_MIXT_VCI2) = "vci2"
  pRegion%plot%pvNameShort(PV_MIXT_VCL2) = "vcl2"
  pRegion%plot%pvNameShort(PV_MIXT_VCLI) = "vcli"
  pRegion%plot%pvNameShort(PV_MIXT_VCLR) = "vclr"

  pRegion%plot%pvNameLong(PV_MIXT_VCI2) = "TBD"  
  pRegion%plot%pvNameLong(PV_MIXT_VCL2) = "TBD"  
  pRegion%plot%pvNameLong(PV_MIXT_VCLI) = "TBD" 
  pRegion%plot%pvNameLong(PV_MIXT_VCLR) = "TBD" 

! ------------------------------------------------------------------------------
! Gradient error 
! ------------------------------------------------------------------------------

  pRegion%plot%pvNameShort(PV_MIXT_GREX) = "grex"
  pRegion%plot%pvNameShort(PV_MIXT_GREY) = "grey"
  pRegion%plot%pvNameShort(PV_MIXT_GREZ) = "grez"

  pRegion%plot%pvNameLong(PV_MIXT_GREX) = "Gradient error x-component"  
  pRegion%plot%pvNameLong(PV_MIXT_GREY) = "Gradient error y-component"  
  pRegion%plot%pvNameLong(PV_MIXT_GREZ) = "Gradient error z-component"         

! ==============================================================================
! Physics modules
! ==============================================================================

#ifdef PLAG
! ------------------------------------------------------------------------------
! Particles 
! ------------------------------------------------------------------------------

  pRegion%plot%pvNameShort(PV_PLAG_DIA3) = "dp3" 
  pRegion%plot%pvNameShort(PV_PLAG_DIA4) = "dp4"
  pRegion%plot%pvNameShort(PV_PLAG_NDNS) = "ndp"
  pRegion%plot%pvNameShort(PV_PLAG_XVEL) = "up"
  pRegion%plot%pvNameShort(PV_PLAG_YVEL) = "vp"
  pRegion%plot%pvNameShort(PV_PLAG_ZVEL) = "wp"
  pRegion%plot%pvNameShort(PV_PLAG_TEMP) = "Tp"
  pRegion%plot%pvNameShort(PV_PLAG_MASS) = "mp"
  
  pRegion%plot%pvNameLong(PV_PLAG_DIA3) = "Particle mean diameter (3)"
  pRegion%plot%pvNameLong(PV_PLAG_DIA4) = "Particle mean diameter (4)"
  pRegion%plot%pvNameLong(PV_PLAG_NDNS) = "Particle number density"
  pRegion%plot%pvNameLong(PV_PLAG_XVEL) = "Particle mean x-velocity"
  pRegion%plot%pvNameLong(PV_PLAG_YVEL) = "Particle mean y-velocity"
  pRegion%plot%pvNameLong(PV_PLAG_ZVEL) = "Particle mean z-velocity"
  pRegion%plot%pvNameLong(PV_PLAG_TEMP) = "Particle mean temperature"
  pRegion%plot%pvNameLong(PV_PLAG_MASS) = "Particle mean mass"                 
#endif

! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BuildPlottingVarNames 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Compute vortex-core plotting variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Assumes that stencils and differentiation weights are accessible.
!
! ******************************************************************************  

SUBROUTINE RFLU_ComputePlottingVarsCore(pRegion)

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

  INTEGER :: errorFlag,icg,iLocVCI2,iLocVCLI,iLocVCLR,iLocVCL2,info,lwork,nDim
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfo
  REAL(RFREAL) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: wi,work,wr
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,o,s,vl,vr
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pPv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ComputePlottingVarsCore', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing vortex-core plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  
  pCv => pRegion%mixt%cv
  pPv => pRegion%plot%pv

  nDim  = pRegion%mixtInput%dimens
  lwork = 3*nDim

! ******************************************************************************
! Get plotting variable locations
! ******************************************************************************

  iLocVCI2 = pRegion%plot%pv2pvi(PV_MIXT_VCI2)
  iLocVCL2 = pRegion%plot%pv2pvi(PV_MIXT_VCL2)
  iLocVCLI = pRegion%plot%pv2pvi(PV_MIXT_VCLI)
  iLocVCLR = pRegion%plot%pv2pvi(PV_MIXT_VCLR)      

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(grad(XCOORD:ZCOORD,3,pGrid%nCellsTot),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
  END IF ! global%error  
  
  ALLOCATE(a(nDim,nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'a')
  END IF ! global%error 
 
  ALLOCATE(o(nDim,nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'o')
  END IF ! global%error     
 
  ALLOCATE(s(nDim,nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'s')
  END IF ! global%error
   
  ALLOCATE(wr(nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'wr')
  END IF ! global%error  
  
  ALLOCATE(wi(nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'wi')
  END IF ! global%error  
  
  ALLOCATE(vl(nDim,nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vl')
  END IF ! global%error  
  
  ALLOCATE(vr(nDim,nDim),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vr')
  END IF ! global%error    

  ALLOCATE(work(lwork),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'work')
  END IF ! global%error 

! ******************************************************************************
! Compute plotting variables
! ******************************************************************************

! ==============================================================================
! Compute velocity gradients
! ==============================================================================

  ALLOCATE(varInfo(CV_MIXT_XVEL:CV_MIXT_ZVEL),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfo')
  END IF ! global%error 

  varInfo(CV_MIXT_XVEL) = CV_MIXT_XVEL
  varInfo(CV_MIXT_YVEL) = CV_MIXT_YVEL
  varInfo(CV_MIXT_ZVEL) = CV_MIXT_ZVEL  

  CALL RFLU_ComputeGradCellsWrapper(pRegion,CV_MIXT_XVEL,CV_MIXT_ZVEL,1,3, &
                                    varInfo,pCv,grad)

  DEALLOCATE(varInfo,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfo')
  END IF ! global%error 
  
! ==============================================================================
! Compute eigenvalues: NOTE does not make sense for 1d.
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%dimens ) 
    CASE ( 2 )   
      DO icg = 1,pGrid%nCellsTot    
        dudx = grad(XCOORD,1,icg)
        dudy = grad(YCOORD,1,icg) 
        
        dvdx = grad(XCOORD,2,icg)
        dvdy = grad(YCOORD,2,icg)               
      
        a(1,1) = dudx
        a(1,2) = dudy
        
        a(2,1) = dvdx
        a(2,2) = dvdy

        CALL dgeev('N','N',nDim,a,nDim,wr,wi,vl,nDim,vr,nDim,work,lwork,info)

        pPv(iLocVCI2,icg) = wr(1)*wr(2)
        pPv(iLocVCLI,icg) = wi(1)
        pPv(iLocVCLR,icg) = wr(1)/(wi(1)+1.0E-12_RFREAL)
        
        s(1,1) = dudx
        s(1,2) = 0.5_RFREAL*(dudy + dvdx)                 
        
        s(2,1) = s(1,2) 
        s(2,2) = dvdy
        
        o(1,1) = 0.0_RFREAL
        o(1,2) = 0.5_RFREAL*(dudy - dvdx)                 
        
        o(2,1) = -o(1,2) 
        o(2,2) = 0.0_RFREAL
        
        a = MATMUL(s,s) + MATMUL(o,o)
                
        CALL dgeev('N','N',nDim,a,nDim,wr,wi,vl,nDim,vr,nDim,work,lwork,info)

        pPv(iLocVCL2,icg) = wr(1)              
      END DO ! icg
    CASE ( 3 ) 
      DO icg = 1,pGrid%nCellsTot
        dudx = grad(XCOORD,1,icg)
        dudy = grad(YCOORD,1,icg)   
        dudz = grad(ZCOORD,1,icg) 

        dvdx = grad(XCOORD,2,icg)
        dvdy = grad(YCOORD,2,icg)   
        dvdz = grad(ZCOORD,2,icg) 

        dwdx = grad(XCOORD,3,icg)
        dwdy = grad(YCOORD,3,icg)   
        dwdz = grad(ZCOORD,3,icg)  
    
        a(1,1) = dudx
        a(1,2) = dudy
        a(1,3) = dudz
        
        a(2,1) = dvdx
        a(2,2) = dvdy
        a(2,3) = dvdz

        a(3,1) = dwdx
        a(3,2) = dwdy
        a(3,3) = dwdz
        
        CALL dgeev('N','N',nDim,a,nDim,wr,wi,vl,nDim,vr,nDim,work,lwork,info)

        pPv(iLocVCI2,icg) = wr(1)*wr(2) + wr(1)*wr(3) + wr(2)*wr(3) 
        pPv(iLocVCLI,icg) = wi(1)
        pPv(iLocVCLR,icg) = wr(1)/(wi(1)+1.0E-12_RFREAL) 
        
        s(1,1) = dudx
        s(1,2) = 0.5_RFREAL*(dudy + dvdx)
        s(1,3) = 0.5_RFREAL*(dudz + dwdx)                         
        
        s(2,1) = s(1,2) 
        s(2,2) = dvdy
        s(2,3) = 0.5_RFREAL*(dvdz + dwdy) 
        
        s(3,1) = s(1,3) 
        s(3,2) = s(2,3)
        s(3,3) = dwdz               
        
        o(1,1) = 0.0_RFREAL
        o(1,2) = 0.5_RFREAL*(dudy - dvdx)
        o(1,3) = 0.5_RFREAL*(dudz - dwdx)                         
        
        o(2,1) = -o(1,2) 
        o(2,2) = 0.0_RFREAL
        o(2,3) = 0.5_RFREAL*(dvdz - dwdy) 
       
        o(3,1) = -o(1,3) 
        o(3,2) = -o(2,3)
        o(3,3) = 0.0_RFREAL                 
        
        a = MATMUL(s,s) + MATMUL(o,o)
                
        CALL dgeev('N','N',nDim,a,nDim,wr,wi,vl,nDim,vr,nDim,work,lwork,info)

        pPv(iLocVCL2,icg) = wr(2)                     
      END DO ! icg    
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%dimens
  
! ==============================================================================
! Deallocate temporary memory
! ==============================================================================

  DEALLOCATE(grad,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
  END IF ! global%error
  
  DEALLOCATE(a,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'a')
  END IF ! global%error  
  
  DEALLOCATE(o,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'o')
  END IF ! global%error  

  DEALLOCATE(s,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'s')
  END IF ! global%error      
  
  DEALLOCATE(wr,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'wr')
  END IF ! global%error  
  
  DEALLOCATE(wi,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'wi')
  END IF ! global%error  
  
  DEALLOCATE(vl,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vl')
  END IF ! global%error  
  
  DEALLOCATE(vr,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vr')
  END IF ! global%error    

  DEALLOCATE(work,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'work')
  END IF ! global%error 

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing vortex-core plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePlottingVarsCore 
  
 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Compute plotting variables related to discontinuities.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Assumes that stencils and differentiation weights are accessible.
!
! ******************************************************************************  

SUBROUTINE RFLU_ComputePlottingVarsDisc(pRegion)

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

  INTEGER :: errorFlag,icg,iLocIntf,iLocSchl,iLocShad,termInt
  INTEGER :: varInfo(CV_MIXT_DENS:CV_MIXT_DENS)
  REAL(RFREAL) :: gradMax,rMax,rMin,term
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pPv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ComputePlottingVarsDisc', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing discontinuity plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  
  pCv => pRegion%mixt%cv
  pPv => pRegion%plot%pv

! ******************************************************************************
! Get and set plotting variable locations 
! ******************************************************************************

  iLocSchl = pRegion%plot%pv2pvi(PV_MIXT_SCHL)
  iLocShad = pRegion%plot%pv2pvi(PV_MIXT_SHAD)
  iLocIntf = pRegion%plot%pv2pvi(PV_MIXT_INTF)    

! ******************************************************************************
! Compute plotting variables used to visualize discontinuities
! ******************************************************************************

! ==============================================================================
! Schlieren and Shadowgraph
! ==============================================================================

! ------------------------------------------------------------------------------
! Allocate temporary memory
! ------------------------------------------------------------------------------

  ALLOCATE(grad(XCOORD:ZCOORD,1,pGrid%nCellsTot),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
  END IF ! global%error  

! ------------------------------------------------------------------------------
! Compute gradients
! ------------------------------------------------------------------------------

  varInfo(CV_MIXT_DENS) = CV_MIXT_DENS  

  CALL RFLU_ComputeGradCellsWrapper(pRegion,CV_MIXT_DENS,CV_MIXT_DENS,1,1, &
                                    varInfo,pCv,grad)

! ------------------------------------------------------------------------------
! Compute Schlieren and Shadowgraph variables. NOTE when using scaled Schlieren,
! will only be computed correctly for single region or merged regions.
! ------------------------------------------------------------------------------

  DO icg = 1,pGrid%nCellsTot
    pPv(iLocSchl,icg) = SQRT(grad(XCOORD,1,icg)**2 & 
                           + grad(YCOORD,1,icg)**2 & 
                           + grad(ZCOORD,1,icg)**2)
    pPv(iLocShad,icg) = 0.0_RFREAL ! TEMPORARY
  END DO ! icg

  IF ( global%postSchType /= 0 ) THEN
    gradMax = MAXVAL(pPv(iLocSchl,1:pGrid%nCellsTot))

    DO icg = 1,pGrid%nCellsTot
      pPv(iLocSchl,icg) = EXP(-global%postSchExp*pPv(iLocSchl,icg)/gradMax)
    END DO ! icg
  END IF ! global%postSchType
  
! ------------------------------------------------------------------------------
! Deallocate temporary memory
! ------------------------------------------------------------------------------

  DEALLOCATE(grad,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
  END IF ! global%error

! ==============================================================================
! Interferogram - NOTE will only be computed correctly for single region or 
! merged regions.
! ==============================================================================

  rMax = MAXVAL(pCv(CV_MIXT_DENS,1:pGrid%nCells))
  rMin = MINVAL(pCv(CV_MIXT_DENS,1:pGrid%nCells))

  DO icg = 1,pGrid%nCellsTot
    term    = (pCv(CV_MIXT_DENS,icg)-rMin)/(rMax-rMin)
    termInt = INT(global%postNFringes*term + 0.5_RFREAL)
  
    pPv(iLocIntf,icg) = MOD(termInt,2) 
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing discontinuity plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePlottingVarsDisc 
 
 
 
 
 

! ******************************************************************************
!
! Purpose: Compute gradient-error plotting variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Assumes that stencils and differentiation weights are accessible.
!
! ******************************************************************************  

SUBROUTINE RFLU_ComputePlottingVarsGradErr(pRegion)

  USE RFLU_ModExactFlow, ONLY: RFLU_SetExactFlowLinear, & 
                               RFLU_SetExactFlowTrig

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

  INTEGER :: errorFlag,icg,iLocXErr,iLocYErr,iLocZErr,iVar,iVarBeg,iVarEnd
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfo
  REAL(RFREAL) :: gxe,gye,gze,nx,ny,nz,var,x,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pPv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ComputePlottingVarsGradErr', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing gradient-error plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  
  pCv => pRegion%mixt%cv
  pPv => pRegion%plot%pv
  
  iVarBeg = CV_MIXT_DENS ! NOTE at present works only for 1 var
  iVarEnd = iVarBeg 

! ******************************************************************************
! Get and set plotting variable locations
! ******************************************************************************

  iLocXErr = pRegion%plot%pv2pvi(PV_MIXT_GREX)
  iLocYErr = pRegion%plot%pv2pvi(PV_MIXT_GREY)
  iLocZErr = pRegion%plot%pv2pvi(PV_MIXT_GREZ)   

! ******************************************************************************
! Compute plotting variables
! ******************************************************************************

! ==============================================================================
! Allocate temporary memory
! ==============================================================================

  ALLOCATE(grad(XCOORD:ZCOORD,iVarBeg:iVarEnd,pGrid%nCellsTot),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
  END IF ! global%error  

! ==============================================================================
! Compute gradients
! ==============================================================================

  ALLOCATE(varInfo(iVarBeg:iVarEnd),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfo')
  END IF ! global%error    

  varInfo(iVarBeg:iVarEnd) = iVarBeg ! NOTE at present works only for 1 var

  CALL RFLU_ComputeGradCellsWrapper(pRegion,iVarBeg,iVarEnd,iVarBeg,iVarEnd, &
                                    varInfo,pCv,grad)

  DEALLOCATE(varInfo,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfo')
  END IF ! global%error    
  
! ==============================================================================
! Compute error
! ==============================================================================

  SELECT CASE ( global%casename ) 
  
! ------------------------------------------------------------------------------
!   Linear function
! ------------------------------------------------------------------------------  
  
    CASE ( "gtlin" )
      SELECT CASE ( pRegion%mixtInput%dimens ) 
        CASE ( 1 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = 0.0_RFREAL
              pPv(iLocZErr,icg) = 0.0_RFREAL   
            END DO ! iVar 
          END DO ! icg
        CASE ( 2 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = grad(YCOORD,iVar,icg)/gye - 1.0_RFREAL
              pPv(iLocZErr,icg) = 0.0_RFREAL   
            END DO ! iVar 
          END DO ! icg
        CASE ( 3 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = grad(YCOORD,iVar,icg)/gye - 1.0_RFREAL
              pPv(iLocZErr,icg) = grad(ZCOORD,iVar,icg)/gze - 1.0_RFREAL  
            END DO ! iVar 
          END DO ! icg        
        CASE DEFAULT ! Defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%dimens           

! ------------------------------------------------------------------------------
!   Trigonometric function
! ------------------------------------------------------------------------------  
  
    CASE ( "gttri" )
      nx = pRegion%mixtInput%prepRealVal1
      ny = pRegion%mixtInput%prepRealVal2
      nz = pRegion%mixtInput%prepRealVal3            
    
      SELECT CASE ( pRegion%mixtInput%dimens ) 
        CASE ( 1 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                         gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = 0.0_RFREAL
              pPv(iLocZErr,icg) = 0.0_RFREAL   
            END DO ! iVar 
          END DO ! icg
        CASE ( 2 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                         gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = grad(YCOORD,iVar,icg)/gye - 1.0_RFREAL
              pPv(iLocZErr,icg) = 0.0_RFREAL   
            END DO ! iVar 
          END DO ! icg
        CASE ( 3 ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                

            DO iVar = iVarBeg,iVarEnd
              CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                         gxe,gye,gze)

              pPv(iLocXErr,icg) = grad(XCOORD,iVar,icg)/gxe - 1.0_RFREAL
              pPv(iLocYErr,icg) = grad(YCOORD,iVar,icg)/gye - 1.0_RFREAL
              pPv(iLocZErr,icg) = grad(ZCOORD,iVar,icg)/gze - 1.0_RFREAL  
            END DO ! iVar 
          END DO ! icg        
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
       END SELECT ! pRegion%mixtInput%dimens           
       
! ------------------------------------------------------------------------------
!   Default
! ------------------------------------------------------------------------------  
         
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%casename
  
! ==============================================================================
! Deallocate temporary memory
! ==============================================================================

  DEALLOCATE(grad,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing gradient-error plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePlottingVarsGradErr 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Compute vorticity plotting variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Assumes that stencils and differentiation weights are accessible.
!
! ******************************************************************************  

SUBROUTINE RFLU_ComputePlottingVarsVort(pRegion)

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

  INTEGER :: errorFlag,icg,iLocXVor,iLocYVor,iLocZVor
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfo
  REAL(RFREAL) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pPv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ComputePlottingVarsVort', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing vorticity plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  
  pCv => pRegion%mixt%cv
  pPv => pRegion%plot%pv

! ******************************************************************************
! Get and set plotting variable locations
! ******************************************************************************

  iLocXVor = pRegion%plot%pv2pvi(PV_MIXT_XVOR)
  iLocYVor = pRegion%plot%pv2pvi(PV_MIXT_YVOR)
  iLocZVor = pRegion%plot%pv2pvi(PV_MIXT_ZVOR)  
        
! ******************************************************************************
! Compute plotting variables
! ******************************************************************************

! ==============================================================================
! Allocate temporary memory
! ==============================================================================

  ALLOCATE(grad(XCOORD:ZCOORD,3,pGrid%nCellsTot),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
  END IF ! global%error  

! ==============================================================================
! Compute gradients
! ==============================================================================

  ALLOCATE(varInfo(CV_MIXT_XVEL:CV_MIXT_ZVEL),STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfo')
  END IF ! global%error 

  varInfo(CV_MIXT_XVEL) = CV_MIXT_XVEL
  varInfo(CV_MIXT_YVEL) = CV_MIXT_YVEL
  varInfo(CV_MIXT_ZVEL) = CV_MIXT_ZVEL  

  CALL RFLU_ComputeGradCellsWrapper(pRegion,CV_MIXT_XVEL,CV_MIXT_ZVEL,1,3, &
                                    varInfo,pCv,grad)

  DEALLOCATE(varInfo,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfo')
  END IF ! global%error 
  
! ==============================================================================
! Compute vorticity
! ==============================================================================

  DO icg = 1,pGrid%nCellsTot
    dudx = grad(XCOORD,1,icg)
    dudy = grad(YCOORD,1,icg)   
    dudz = grad(ZCOORD,1,icg) 
    
    dvdx = grad(XCOORD,2,icg)
    dvdy = grad(YCOORD,2,icg)   
    dvdz = grad(ZCOORD,2,icg) 
    
    dwdx = grad(XCOORD,3,icg)
    dwdy = grad(YCOORD,3,icg)   
    dwdz = grad(ZCOORD,3,icg)         
  
    pPv(iLocXVor,icg) = dwdy - dvdz
    pPv(iLocYVor,icg) = dudz - dwdx
    pPv(iLocZVor,icg) = dvdx - dudy    
  END DO ! icg

! ==============================================================================
! Deallocate temporary memory
! ==============================================================================

  DEALLOCATE(grad,STAT=errorFlag) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing vorticity plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePlottingVarsVort 
 
 
 
 
 
 



! ******************************************************************************
!
! Purpose: Wrapper for computing plotting variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. Assumes that stencils and differentiation weights are accessible.
!
! ******************************************************************************  

SUBROUTINE RFLU_ComputePlottingVarsWrapper(pRegion)

  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, & 
                               RFLU_ConvertCvPrim2Cons

#ifdef PLAG
  USE PLAG_ModEulerian
#endif

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

  INTEGER :: icgBeg,icgEnd,iPvBeg,iPvEnd
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ComputePlottingVarsWrapper', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Mixture
! ******************************************************************************

  IF ( global%postDiscFlag .EQV. .TRUE. ) THEN   
    CALL RFLU_ComputePlottingVarsDisc(pRegion)
  END IF ! global%postDiscFlag

  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

  IF ( global%postGradFlag .EQV. .TRUE. ) THEN   
    CALL RFLU_ComputePlottingVarsGradErr(pRegion)
  END IF ! global%postGradFlag 

  IF ( global%postVortFlag .EQV. .TRUE. ) THEN   
    CALL RFLU_ComputePlottingVarsVort(pRegion)
  END IF ! global%postVortFlag 
  
  IF ( global%postVortCoreFlag .EQV. .TRUE. ) THEN   
    CALL RFLU_ComputePlottingVarsCore(pRegion)
  END IF ! global%postVortCoreFlag   

  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

! ******************************************************************************
! Physics modules
! ******************************************************************************

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    IF ( global%postLag2EulFlag .EQV. .TRUE. ) THEN
      iPvBeg = pRegion%plot%pv2pvi(PV_PLAG_DIA3)
      iPvEnd = pRegion%plot%pv2pvi(PV_PLAG_MASS)
      icgBeg = 1
      icgEnd = pRegion%grid%nCells
  
      CALL PLAG_RFLU_CalcEulerianField(pRegion, & 
           pRegion%plot%pv(iPvBeg:iPvEnd,icgBeg:icgEnd))
    END IF ! global%postLag2EulFlag
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputePlottingVarsWrapper
 
 
 
 
 
 
 
 
 
 
! ******************************************************************************
!
! Purpose: Count plotting variables.
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

SUBROUTINE RFLU_CountPlottingVars(pRegion)

#ifdef PLAG
  USE PLAG_ModParameters
#endif

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

  INTEGER :: errorFlag,iPv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_CountPlottingVars', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Initialize 
! ******************************************************************************

  DO iPv = 1,PV_XXXX_NVAR
    pRegion%plot%pv2pvi(iPv) = CRAZY_VALUE_INT ! Important
  END DO ! iPv

! ******************************************************************************
! Count number of plotting variables
! ******************************************************************************

  pRegion%plot%nPv = 0

! ==============================================================================
! Mixture
! ==============================================================================
        
! ------------------------------------------------------------------------------
! Variables to visualize discontinuities
! ------------------------------------------------------------------------------
  
  IF ( global%postDiscFlag .EQV. .TRUE. ) THEN
    pRegion%plot%pv2pvi(PV_MIXT_SCHL) = pRegion%plot%nPv + 1
    pRegion%plot%pv2pvi(PV_MIXT_SHAD) = pRegion%plot%nPv + 2
    pRegion%plot%pv2pvi(PV_MIXT_INTF) = pRegion%plot%nPv + 3     
  
    pRegion%plot%nPv = pRegion%plot%nPv + 3
  END IF ! global%postDiscFlag

! ------------------------------------------------------------------------------
! Vorticity
! ------------------------------------------------------------------------------

  IF ( global%postVortFlag .EQV. .TRUE. ) THEN
    pRegion%plot%pv2pvi(PV_MIXT_XVOR) = pRegion%plot%nPv + 1
    pRegion%plot%pv2pvi(PV_MIXT_YVOR) = pRegion%plot%nPv + 2
    pRegion%plot%pv2pvi(PV_MIXT_ZVOR) = pRegion%plot%nPv + 3     
  
    pRegion%plot%nPv = pRegion%plot%nPv + 3
  END IF ! global%postVortFlag
  
! ------------------------------------------------------------------------------
! Vortex core
! ------------------------------------------------------------------------------

  IF ( global%postVortCoreFlag .EQV. .TRUE. ) THEN
    pRegion%plot%pv2pvi(PV_MIXT_VCI2) = pRegion%plot%nPv + 1
    pRegion%plot%pv2pvi(PV_MIXT_VCL2) = pRegion%plot%nPv + 2
    pRegion%plot%pv2pvi(PV_MIXT_VCLI) = pRegion%plot%nPv + 3
    pRegion%plot%pv2pvi(PV_MIXT_VCLR) = pRegion%plot%nPv + 4                 
  
    pRegion%plot%nPv = pRegion%plot%nPv + 4
  END IF ! global%postVortCoreFlag  

! ------------------------------------------------------------------------------
! Variables to visualize gradient error
! ------------------------------------------------------------------------------
  
  IF ( global%postGradFlag .EQV. .TRUE. ) THEN
    pRegion%plot%pv2pvi(PV_MIXT_GREX) = pRegion%plot%nPv + 1
    pRegion%plot%pv2pvi(PV_MIXT_GREY) = pRegion%plot%nPv + 2
    pRegion%plot%pv2pvi(PV_MIXT_GREZ) = pRegion%plot%nPv + 3     
  
    pRegion%plot%nPv = pRegion%plot%nPv + 3
  END IF ! global%postGradFlag    

! ==============================================================================
! Physics modules
! ==============================================================================

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN 
    IF ( global%postLag2EulFlag .EQV. .TRUE. ) THEN
      pRegion%plot%pv2pvi(PV_PLAG_DIA3) = pRegion%plot%nPv + 1 
      pRegion%plot%pv2pvi(PV_PLAG_DIA4) = pRegion%plot%nPv + 2 
      pRegion%plot%pv2pvi(PV_PLAG_NDNS) = pRegion%plot%nPv + 3
      pRegion%plot%pv2pvi(PV_PLAG_XVEL) = pRegion%plot%nPv + 4
      pRegion%plot%pv2pvi(PV_PLAG_YVEL) = pRegion%plot%nPv + 5
      pRegion%plot%pv2pvi(PV_PLAG_ZVEL) = pRegion%plot%nPv + 6
      pRegion%plot%pv2pvi(PV_PLAG_TEMP) = pRegion%plot%nPv + 7
      pRegion%plot%pv2pvi(PV_PLAG_MASS) = pRegion%plot%nPv + 8               
  
      pRegion%plot%nPv = pRegion%plot%nPv + EV_PLAG_LAST
    END IF ! global%postLag2EulFlag
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! Print information and build names
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME, &
                                   'Number of plotting variables:', &
                                   pRegion%plot%nPv
  END IF ! global%verbLevel

  CALL RFLU_BuildPlottingVarNames(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Counting plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CountPlottingVars
  
  
  
  
 
 
 
! ******************************************************************************
!
! Purpose: Create plotting variable maps.
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

SUBROUTINE RFLU_CreatePlottingVarMaps(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CreatePlottingVarMaps', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating plotting variable maps...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    ALLOCATE(pRegion%plot%pvi2pv(pRegion%plot%nPv),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%plot%pvi2pv')
    END IF ! global%error  
  ELSE
    CALL RFLU_NullifyPlottingVarMaps(pRegion)
  END IF ! pRegion%plot%nPv

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating plotting variable maps done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreatePlottingVarMaps 
 
 
  
  
  


! ******************************************************************************
!
! Purpose: Create plotting variables.
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

SUBROUTINE RFLU_CreatePlottingVars(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CreatePlottingVars', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    ALLOCATE(pRegion%plot%pv(pRegion%plot%nPv,pRegion%grid%nCellsTot), & 
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%plot%pv')
    END IF ! global%error

    CALL RFLU_InitPlottingVars(pRegion)
  ELSE
    CALL RFLU_NullifyPlottingVars(pRegion)
  END IF ! pRegion%plot%nPv

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreatePlottingVars  
  
  
  





! ******************************************************************************
!
! Purpose: Create plotting variables at vertices.
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

SUBROUTINE RFLU_CreatePlottingVarsVert(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CreatePlottingVarsVert', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Create plotting variables
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    ALLOCATE(pRegion%plot%pvVert(pRegion%plot%nPv,pRegion%grid%nVertTot), & 
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%plot%pvVert')
    END IF ! global%error

    CALL RFLU_InitPlottingVarsVert(pRegion)
  ELSE
    CALL RFLU_NullifyPlottingVarsVert(pRegion)
  END IF ! pRegion%plot%nPv

! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreatePlottingVarsVert








! ******************************************************************************
!
! Purpose: Decide whether plotting vars to be computed.
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

LOGICAL FUNCTION RFLU_DecideComputePlottingVars(pRegion)

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

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

! ******************************************************************************
! Initialize
! ******************************************************************************

  RFLU_DecideComputePlottingVars = .FALSE.

! ******************************************************************************
! Decide whether plotting vars to be computed
! ******************************************************************************

  IF ( pRegion%mixtInput%fluidModel == FLUID_MODEL_COMP ) THEN
    IF ( (global%postDiscFlag     .EQV. .TRUE.) .OR. &
         (global%postGradFlag     .EQV. .TRUE.) .OR. & 
         (global%postVortFlag     .EQV. .TRUE.) .OR. &
         (global%postVortCoreFlag .EQV. .TRUE.) .OR. & 
         (global%postLag2EulFlag  .EQV. .TRUE.) ) THEN
      RFLU_DecideComputePlottingVars = .TRUE.
    END IF ! global%postDiscFlag
  END IF ! pRegion%mixtInput%fluidModel
  
! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_DecideComputePlottingVars 







! ******************************************************************************
!
! Purpose: Destroy plotting variable maps.
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

SUBROUTINE RFLU_DestroyPlottingVarMaps(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyPlottingVarMaps', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying plotting variable maps...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    DEALLOCATE(pRegion%plot%pvi2pv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%plot%pvi2pv')
    END IF ! global%error  
  END IF ! pRegion%plot%nPv

  CALL RFLU_NullifyPlottingVarMaps(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Destroying plotting variable maps done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyPlottingVarMaps 







! ******************************************************************************
!
! Purpose: Destroy plotting variables.
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

SUBROUTINE RFLU_DestroyPlottingVars(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyPlottingVars', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying plotting variables...'
  END IF ! global%verbLevel
  
! ******************************************************************************
! Deallocate memory
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    DEALLOCATE(pRegion%plot%pv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%plot%pv')
    END IF ! global%error
  END IF ! pRegion%plot%nPv
  
  CALL RFLU_NullifyPlottingVars(pRegion)
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyPlottingVars 








! ******************************************************************************
!
! Purpose: Destroy plotting variables at vertices.
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

SUBROUTINE RFLU_DestroyPlottingVarsVert(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DestroyPlottingVarsVert', & 
                        'RFLU_ModPlottingVars.F90')
  
! ******************************************************************************
! Destroy plotting variables
! ******************************************************************************
  
  IF ( pRegion%plot%nPv > 0 ) THEN 
    DEALLOCATE(pRegion%plot%pvVert,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%plot%pvVert')
    END IF ! global%error
  END IF ! pRegion%plot%nPv
  
  CALL RFLU_NullifyPlottingVarsVert(pRegion)
  
! ******************************************************************************
! End
! ******************************************************************************
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DestroyPlottingVarsVert 






! ******************************************************************************
!
! Purpose: Initialize plotting variables.
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

SUBROUTINE RFLU_InitPlottingVars(pRegion)

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

  INTEGER :: icg,iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitPlottingVars', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Initialize plotting variables
! ******************************************************************************
  
  DO iVar = 1,pRegion%plot%nPv
    DO icg = 1,pRegion%grid%nCellsTot
      pRegion%plot%pv(iVar,icg) = 0.0_RFREAL
    END DO ! icg
  END DO ! iVar

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitPlottingVars  







! ******************************************************************************
!
! Purpose: Initialize vertex plotting variables.
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

SUBROUTINE RFLU_InitPlottingVarsVert(pRegion)

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

  INTEGER :: ivg,iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitPlottingVarsVert', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Initializing vertex plotting variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Initialize vertex plotting variables
! ******************************************************************************
  
  DO iVar = 1,pRegion%plot%nPv
    DO ivg = 1,pRegion%grid%nVertTot
      pRegion%plot%pvVert(iVar,ivg) = 0.0_RFREAL
    END DO ! ivg
  END DO ! iVar

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing vertex plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitPlottingVarsVert







! ******************************************************************************
!
! Purpose: Nullify plotting variable maps.
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

SUBROUTINE RFLU_NullifyPlottingVarMaps(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NullifyPlottingVarMaps', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Nullify plotting variables
! ******************************************************************************
  
  NULLIFY(pRegion%plot%pvi2pv)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NullifyPlottingVarMaps





! ******************************************************************************
!
! Purpose: Nullify plotting variables.
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

SUBROUTINE RFLU_NullifyPlottingVars(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NullifyPlottingVars', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Nullify plotting variables
! ******************************************************************************
  
  NULLIFY(pRegion%plot%pv)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NullifyPlottingVars






! ******************************************************************************
!
! Purpose: Nullify plotting variables at vertices.
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

SUBROUTINE RFLU_NullifyPlottingVarsVert(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NullifyPlottingVarsVert', & 
                        'RFLU_ModPlottingVars.F90')

! ******************************************************************************
! Nullify plotting variables
! ******************************************************************************
  
  NULLIFY(pRegion%plot%pvVert)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NullifyPlottingVarsVert






! ******************************************************************************
!
! Purpose: Print plotting variable info.
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

SUBROUTINE RFLU_PrintPlottingVarsInfo(pRegion)

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

  INTEGER :: errorFlag,iPv,iPv2
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_PrintPlottingVarsInfo', & 
                        'RFLU_ModPlottingVars.F90')

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Printing information on plotting variables...'
  END IF ! global%verbLevel
  
! ******************************************************************************
! Print names of plotting variables
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,5X,A,2X,A,2X,A,26X,A)') SOLVER_NAME,'#','Index', & 
             'Long name','Short name'

    DO iPv = 1,pRegion%plot%nPv
      iPv2 = pRegion%plot%pvi2pv(iPv)
    
      WRITE(STDOUT,'(A,4X,I2,4X,I2,3X,A32,3X,A6)') SOLVER_NAME,iPv,iPv2, &
            ADJUSTL(pRegion%plot%pvNameLong(iPv2)), & 
            ADJUSTL(pRegion%plot%pvNameShort(iPv2))
    END DO ! iPv
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************
  
  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Printing information on plotting variables done.'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintPlottingVarsInfo





END MODULE RFLU_ModPlottingVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPlottingVars.F90,v $
! Revision 1.15  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2007/03/19 21:44:05  haselbac
! Significant reorganization with move of pv from mixt to plot
!
! Revision 1.12  2007/02/27 13:24:36  haselbac
! Enabled 1d computations
!
! Revision 1.11  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.10  2006/01/06 22:20:15  haselbac
! Added routine to compute grad error
!
! Revision 1.9  2005/12/10 16:56:36  haselbac
! Added decide and init routines
!
! Revision 1.8  2005/11/30 22:25:25  fnajjar
! Added call to compute Eulerian vars for PLAG
!
! Revision 1.7  2005/10/05 16:20:49  haselbac
! Bug fix in USE statement
!
! Revision 1.6  2005/08/10 00:38:18  haselbac
! Added routine for vortex core extraction and wrapper
!
! Revision 1.5  2005/07/25 12:24:32  haselbac
! Added vorticity plotting vars, moved stencil and weights out for efficiency 
! and commonality
!
! Revision 1.4  2005/07/05 19:24:03  haselbac
! Added exponential Schlieren indicator
!
! Revision 1.3  2005/07/01 21:29:28  haselbac
! Bug fix: Allocate pv memory for nCellsTot
!
! Revision 1.2  2005/05/03 20:38:22  haselbac
! Removed DEBUG code
!
! Revision 1.1  2005/05/01 14:17:11  haselbac
! Initial revision
!
! ******************************************************************************


























