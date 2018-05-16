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
! Purpose: Compute grid speeds.
!
! Description: The grid speeds are computed by computing the volume swept out
!   by the motion of the face, divided by a surface area and the time step. 
!   Depending on whether the face is triangular or quadrilateral, the volume
!   swept out by the face motion is a prism or a hexahedron.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. The order in which RFLU_ComputeGridSpeeds and RFLU_BuildGeometry are 
!      called is important because the grid speeds are defined to be a volume
!      change divided by a face area (and the time-step). If the face area is
!      from the old grid, the grid speeds will not satisfy continuity when 
!      used to compute fluxes.
!   2. The volumes of the polyhedra swept out by the triangular and 
!      quadrilateral faces are computed in a manner consistent with the way in
!      which the cell volumes are computed, see, e.g., RFLU_ModGeometry.F90. 
!
! *****************************************************************************
!
! $Id: RFLU_ComputeGridSpeeds.F90,v 1.9 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeGridSpeeds(pRegion)

  USE ModDataTypes
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  USE RFLU_ModGrid

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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,ifc,ifcl,iPatch,iv,v1,v1l,v2,v2l,v3,v3l,v4,v4l
  REAL(RFREAL), PARAMETER :: THRD = 1.0_RFREAL/3.0_RFREAL  
  REAL(RFREAL) :: fcx,fcy,fcz,fnx,fny,fnz,term
  REAL(RFREAL) :: xyzAvg(XCOORD:ZCOORD)  
  REAL(RFREAL) :: xyzHex(XCOORD:ZCOORD,8),xyzNodes(XCOORD:ZCOORD,4), & 
                  xyzPri(XCOORD:ZCOORD,6)
  REAL(RFREAL) :: vsnew, vsold
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pXyz,pXyzOld,pXyzOld2
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld,pGridOld2
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeGridSpeeds.F90,v $ $Revision: 1.9 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeGridSpeeds',&
  'RFLU_ComputeGridSpeeds.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid speeds...'
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid speed extrema:'    
  END IF ! global%myProcid 
  
! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
    
  pGrid     => pRegion%grid
  pGridOld  => pRegion%gridOld
  pGridOld2 => pRegion%gridOld2
  
  pXyz      => pGrid%xyz
  pXyzOld   => pGridOld%xyz
  pXyzOld2  => pGridOld2%xyz

! ******************************************************************************
! Interior faces
! ******************************************************************************

  DO ifc = 1,pGrid%nFaces 
    v1 = pGrid%f2v(1,ifc)
    v2 = pGrid%f2v(2,ifc)
    v3 = pGrid%f2v(3,ifc)
    v4 = pGrid%f2v(4,ifc)
             
    pGrid%gs(ifc) = 0.0_RFREAL

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
       vsnew = 0.0_RFREAL
       vsold = 0.0_RFREAL
    END IF ! global%solverType
      
! ==============================================================================
!   Triangular face, sweeps out prism during grid motion
! ==============================================================================    
    
    IF ( v4 == VERT_NONE ) THEN ! triangular face
      xyzPri(XCOORD,1) = pXyzOld(XCOORD,v1)
      xyzPri(XCOORD,2) = pXyzOld(XCOORD,v2)
      xyzPri(XCOORD,3) = pXyzOld(XCOORD,v3)            
    
      xyzPri(YCOORD,1) = pXyzOld(YCOORD,v1)
      xyzPri(YCOORD,2) = pXyzOld(YCOORD,v2)
      xyzPri(YCOORD,3) = pXyzOld(YCOORD,v3)     
    
      xyzPri(ZCOORD,1) = pXyzOld(ZCOORD,v1)
      xyzPri(ZCOORD,2) = pXyzOld(ZCOORD,v2)
      xyzPri(ZCOORD,3) = pXyzOld(ZCOORD,v3)      
    
      xyzPri(XCOORD,4) = pXyz(XCOORD,v1)
      xyzPri(XCOORD,5) = pXyz(XCOORD,v2)
      xyzPri(XCOORD,6) = pXyz(XCOORD,v3)            
    
      xyzPri(YCOORD,4) = pXyz(YCOORD,v1)
      xyzPri(YCOORD,5) = pXyz(YCOORD,v2)
      xyzPri(YCOORD,6) = pXyz(YCOORD,v3)     
    
      xyzPri(ZCOORD,4) = pXyz(ZCOORD,v1)
      xyzPri(ZCOORD,5) = pXyz(ZCOORD,v2)
      xyzPri(ZCOORD,6) = pXyz(ZCOORD,v3)      
          
! ------------------------------------------------------------------------------
!     Compute average coordinates (approximate centroid)
! ------------------------------------------------------------------------------

      xyzAvg(XCOORD) = 0.0_RFREAL
      xyzAvg(YCOORD) = 0.0_RFREAL
      xyzAvg(ZCOORD) = 0.0_RFREAL            

      term = 1.0_RFREAL/6.0_RFREAL

      DO iv = 1,6
        xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzPri(XCOORD,iv)
        xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzPri(YCOORD,iv)
        xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzPri(ZCOORD,iv)                
      END DO ! iv   
          
      xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
      xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
      xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                      
          
! ------------------------------------------------------------------------------
!     Loop over the five faces of the prism to compute volume
! ------------------------------------------------------------------------------    
    
      DO ifcl = 1,5
        v1l = f2vPri(1,ifcl)
        v2l = f2vPri(2,ifcl)
        v3l = f2vPri(3,ifcl)
        v4l = f2vPri(4,ifcl)        
        
        IF ( v4l == 0 ) THEN ! triangular face
          xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
          xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
          xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)

          xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
          xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
          xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)

          xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
          xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
          xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)

          CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fnx,fny,fnz)
          CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fcx,fcy,fcz)        
        ELSE 
          xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
          xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
          xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
          xyzNodes(XCOORD,4) = xyzPri(XCOORD,v4l) - xyzAvg(XCOORD)          

          xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
          xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
          xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
          xyzNodes(YCOORD,4) = xyzPri(YCOORD,v4l) - xyzAvg(YCOORD)

          xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
          xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
          xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
          xyzNodes(ZCOORD,4) = xyzPri(ZCOORD,v4l) - xyzAvg(ZCOORD)                
        
          CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
          CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
        END IF ! v4l                           
                       
        IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
           vsnew = vsnew + (fcx*fnx + fcy*fny + fcz*fnz)
        ELSE
           pGrid%gs(ifc) = pGrid%gs(ifc) + (fcx*fnx + fcy*fny + fcz*fnz)          
        END IF ! global%solverType

      END DO ! ifcl
    
! ==============================================================================
!   Quadrilateral face, sweeps out hexahedron during grid motion
! ==============================================================================    
    
    ELSE ! quadrilateral face
      xyzHex(XCOORD,1) = pXyzOld(XCOORD,v1)
      xyzHex(XCOORD,2) = pXyzOld(XCOORD,v2)
      xyzHex(XCOORD,3) = pXyzOld(XCOORD,v3)            
      xyzHex(XCOORD,4) = pXyzOld(XCOORD,v4)    
    
      xyzHex(YCOORD,1) = pXyzOld(YCOORD,v1)
      xyzHex(YCOORD,2) = pXyzOld(YCOORD,v2)
      xyzHex(YCOORD,3) = pXyzOld(YCOORD,v3)     
      xyzHex(YCOORD,4) = pXyzOld(YCOORD,v4)  
    
      xyzHex(ZCOORD,1) = pXyzOld(ZCOORD,v1)
      xyzHex(ZCOORD,2) = pXyzOld(ZCOORD,v2)
      xyzHex(ZCOORD,3) = pXyzOld(ZCOORD,v3)
      xyzHex(ZCOORD,4) = pXyzOld(ZCOORD,v4)              
    
      xyzHex(XCOORD,5) = pXyz(XCOORD,v1)
      xyzHex(XCOORD,6) = pXyz(XCOORD,v2)
      xyzHex(XCOORD,7) = pXyz(XCOORD,v3) 
      xyzHex(XCOORD,8) = pXyz(XCOORD,v4)                  
    
      xyzHex(YCOORD,5) = pXyz(YCOORD,v1)
      xyzHex(YCOORD,6) = pXyz(YCOORD,v2)
      xyzHex(YCOORD,7) = pXyz(YCOORD,v3)     
      xyzHex(YCOORD,8) = pXyz(YCOORD,v4)   
    
      xyzHex(ZCOORD,5) = pXyz(ZCOORD,v1)
      xyzHex(ZCOORD,6) = pXyz(ZCOORD,v2)
      xyzHex(ZCOORD,7) = pXyz(ZCOORD,v3)   
      xyzHex(ZCOORD,8) = pXyz(ZCOORD,v4)            
    
! ------------------------------------------------------------------------------
!     Compute average coordinates (approximate centroid)
! ------------------------------------------------------------------------------

      xyzAvg(XCOORD) = 0.0_RFREAL
      xyzAvg(YCOORD) = 0.0_RFREAL
      xyzAvg(ZCOORD) = 0.0_RFREAL            

      term = 1.0_RFREAL/8.0_RFREAL

      DO iv = 1,8
        xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzHex(XCOORD,iv)
        xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzHex(YCOORD,iv)
        xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzHex(ZCOORD,iv)                
      END DO ! iv       
    
      xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
      xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
      xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                          
    
! ------------------------------------------------------------------------------
!     Loop over the six faces of the prism to compute volume
! ------------------------------------------------------------------------------    
        
      DO ifcl = 1,6
        v1l = f2vHex(1,ifcl)
        v2l = f2vHex(2,ifcl)
        v3l = f2vHex(3,ifcl)
        v4l = f2vHex(4,ifcl)
        
        xyzNodes(XCOORD,1) = xyzHex(XCOORD,v1l) - xyzAvg(XCOORD)
        xyzNodes(XCOORD,2) = xyzHex(XCOORD,v2l) - xyzAvg(XCOORD)
        xyzNodes(XCOORD,3) = xyzHex(XCOORD,v3l) - xyzAvg(XCOORD)
        xyzNodes(XCOORD,4) = xyzHex(XCOORD,v4l) - xyzAvg(XCOORD)       

        xyzNodes(YCOORD,1) = xyzHex(YCOORD,v1l) - xyzAvg(YCOORD)
        xyzNodes(YCOORD,2) = xyzHex(YCOORD,v2l) - xyzAvg(YCOORD)
        xyzNodes(YCOORD,3) = xyzHex(YCOORD,v3l) - xyzAvg(YCOORD)
        xyzNodes(YCOORD,4) = xyzHex(YCOORD,v4l) - xyzAvg(YCOORD)       

        xyzNodes(ZCOORD,1) = xyzHex(ZCOORD,v1l) - xyzAvg(ZCOORD)
        xyzNodes(ZCOORD,2) = xyzHex(ZCOORD,v2l) - xyzAvg(ZCOORD)
        xyzNodes(ZCOORD,3) = xyzHex(ZCOORD,v3l) - xyzAvg(ZCOORD)
        xyzNodes(ZCOORD,4) = xyzHex(ZCOORD,v4l) - xyzAvg(ZCOORD)       

        CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
        CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
            
        IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
           vsnew = vsnew + (fcx*fnx + fcy*fny + fcz*fnz)
        ELSE
           pGrid%gs(ifc) = pGrid%gs(ifc) + (fcx*fnx + fcy*fny + fcz*fnz)
        END IF ! global%solverType

      END DO ! ifcl      
    END IF ! v4

! ******************************************************************************
! Interior faces with Old2 grid for implicit code
! ******************************************************************************

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN

! ==============================================================================
!   Triangular face, sweeps out prism during grid motion
! ==============================================================================    
    
       IF ( v4 == VERT_NONE ) THEN ! triangular face
          xyzPri(XCOORD,1) = pXyzOld2(XCOORD,v1)
          xyzPri(XCOORD,2) = pXyzOld2(XCOORD,v2)
          xyzPri(XCOORD,3) = pXyzOld2(XCOORD,v3)            
          
          xyzPri(YCOORD,1) = pXyzOld2(YCOORD,v1)
          xyzPri(YCOORD,2) = pXyzOld2(YCOORD,v2)
          xyzPri(YCOORD,3) = pXyzOld2(YCOORD,v3)     
          
          xyzPri(ZCOORD,1) = pXyzOld2(ZCOORD,v1)
          xyzPri(ZCOORD,2) = pXyzOld2(ZCOORD,v2)
          xyzPri(ZCOORD,3) = pXyzOld2(ZCOORD,v3)      
          
          xyzPri(XCOORD,4) = pXyzOld(XCOORD,v1)
          xyzPri(XCOORD,5) = pXyzOld(XCOORD,v2)
          xyzPri(XCOORD,6) = pXyzOld(XCOORD,v3)            
          
          xyzPri(YCOORD,4) = pXyzOld(YCOORD,v1)
          xyzPri(YCOORD,5) = pXyzOld(YCOORD,v2)
          xyzPri(YCOORD,6) = pXyzOld(YCOORD,v3)     
          
          xyzPri(ZCOORD,4) = pXyzOld(ZCOORD,v1)
          xyzPri(ZCOORD,5) = pXyzOld(ZCOORD,v2)
          xyzPri(ZCOORD,6) = pXyzOld(ZCOORD,v3)      
          
! ------------------------------------------------------------------------------
!     Compute average coordinates (approximate centroid)
! ------------------------------------------------------------------------------

          xyzAvg(XCOORD) = 0.0_RFREAL
          xyzAvg(YCOORD) = 0.0_RFREAL
          xyzAvg(ZCOORD) = 0.0_RFREAL            
          
          term = 1.0_RFREAL/6.0_RFREAL
          
          DO iv = 1,6
             xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzPri(XCOORD,iv)
             xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzPri(YCOORD,iv)
             xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzPri(ZCOORD,iv)                
          END DO ! iv   
          
          xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
          xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
          xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                      
          
! ------------------------------------------------------------------------------
!     Loop over the five faces of the prism to compute volume
! ------------------------------------------------------------------------------    
    
          DO ifcl = 1,5
             v1l = f2vPri(1,ifcl)
             v2l = f2vPri(2,ifcl)
             v3l = f2vPri(3,ifcl)
             v4l = f2vPri(4,ifcl)        
             
             IF ( v4l == 0 ) THEN ! triangular face
                xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
                xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
                xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
                
                xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
                xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
                xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
                
                xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
                xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
                xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
                
                CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fnx,fny,fnz)
                CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fcx,fcy,fcz)        
             ELSE 
                xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
                xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
                xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
                xyzNodes(XCOORD,4) = xyzPri(XCOORD,v4l) - xyzAvg(XCOORD)          
                
                xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
                xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
                xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
                xyzNodes(YCOORD,4) = xyzPri(YCOORD,v4l) - xyzAvg(YCOORD)
                
                xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
                xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
                xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
                xyzNodes(ZCOORD,4) = xyzPri(ZCOORD,v4l) - xyzAvg(ZCOORD)                
                
                CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
                CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
             END IF ! v4l                           
             
             vsold = vsold + (fcx*fnx + fcy*fny + fcz*fnz)
             
          END DO ! ifcl
    
! ==============================================================================
!   Quadrilateral face, sweeps out hexahedron during grid motion
! ==============================================================================    
    
       ELSE ! quadrilateral face
          xyzHex(XCOORD,1) = pXyzOld2(XCOORD,v1)
          xyzHex(XCOORD,2) = pXyzOld2(XCOORD,v2)
          xyzHex(XCOORD,3) = pXyzOld2(XCOORD,v3)            
          xyzHex(XCOORD,4) = pXyzOld2(XCOORD,v4)    
          
          xyzHex(YCOORD,1) = pXyzOld2(YCOORD,v1)
          xyzHex(YCOORD,2) = pXyzOld2(YCOORD,v2)
          xyzHex(YCOORD,3) = pXyzOld2(YCOORD,v3)     
          xyzHex(YCOORD,4) = pXyzOld2(YCOORD,v4)  
          
          xyzHex(ZCOORD,1) = pXyzOld2(ZCOORD,v1)
          xyzHex(ZCOORD,2) = pXyzOld2(ZCOORD,v2)
          xyzHex(ZCOORD,3) = pXyzOld2(ZCOORD,v3)
          xyzHex(ZCOORD,4) = pXyzOld2(ZCOORD,v4)              
          
          xyzHex(XCOORD,5) = pXyzOld(XCOORD,v1)
          xyzHex(XCOORD,6) = pXyzOld(XCOORD,v2)
          xyzHex(XCOORD,7) = pXyzOld(XCOORD,v3) 
          xyzHex(XCOORD,8) = pXyzOld(XCOORD,v4)                  
          
          xyzHex(YCOORD,5) = pXyzOld(YCOORD,v1)
          xyzHex(YCOORD,6) = pXyzOld(YCOORD,v2)
          xyzHex(YCOORD,7) = pXyzOld(YCOORD,v3)     
          xyzHex(YCOORD,8) = pXyzOld(YCOORD,v4)   
          
          xyzHex(ZCOORD,5) = pXyzOld(ZCOORD,v1)
          xyzHex(ZCOORD,6) = pXyzOld(ZCOORD,v2)
          xyzHex(ZCOORD,7) = pXyzOld(ZCOORD,v3)   
          xyzHex(ZCOORD,8) = pXyzOld(ZCOORD,v4)            
          
! ------------------------------------------------------------------------------
!     Compute average coordinates (approximate centroid)
! ------------------------------------------------------------------------------

          xyzAvg(XCOORD) = 0.0_RFREAL
          xyzAvg(YCOORD) = 0.0_RFREAL
          xyzAvg(ZCOORD) = 0.0_RFREAL            
          
          term = 1.0_RFREAL/8.0_RFREAL
          
          DO iv = 1,8
             xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzHex(XCOORD,iv)
             xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzHex(YCOORD,iv)
             xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzHex(ZCOORD,iv)                
          END DO ! iv       
          
          xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
          xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
          xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                          
    
! ------------------------------------------------------------------------------
!     Loop over the six faces of the prism to compute volume
! ------------------------------------------------------------------------------    
        
          DO ifcl = 1,6
             v1l = f2vHex(1,ifcl)
             v2l = f2vHex(2,ifcl)
             v3l = f2vHex(3,ifcl)
             v4l = f2vHex(4,ifcl)
             
             xyzNodes(XCOORD,1) = xyzHex(XCOORD,v1l) - xyzAvg(XCOORD)
             xyzNodes(XCOORD,2) = xyzHex(XCOORD,v2l) - xyzAvg(XCOORD)
             xyzNodes(XCOORD,3) = xyzHex(XCOORD,v3l) - xyzAvg(XCOORD)
             xyzNodes(XCOORD,4) = xyzHex(XCOORD,v4l) - xyzAvg(XCOORD)       
             
             xyzNodes(YCOORD,1) = xyzHex(YCOORD,v1l) - xyzAvg(YCOORD)
             xyzNodes(YCOORD,2) = xyzHex(YCOORD,v2l) - xyzAvg(YCOORD)
             xyzNodes(YCOORD,3) = xyzHex(YCOORD,v3l) - xyzAvg(YCOORD)
             xyzNodes(YCOORD,4) = xyzHex(YCOORD,v4l) - xyzAvg(YCOORD)       
             
             xyzNodes(ZCOORD,1) = xyzHex(ZCOORD,v1l) - xyzAvg(ZCOORD)
             xyzNodes(ZCOORD,2) = xyzHex(ZCOORD,v2l) - xyzAvg(ZCOORD)
             xyzNodes(ZCOORD,3) = xyzHex(ZCOORD,v3l) - xyzAvg(ZCOORD)
             xyzNodes(ZCOORD,4) = xyzHex(ZCOORD,v4l) - xyzAvg(ZCOORD)       
             
             CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
             CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
             
             vsold = vsold + (fcx*fnx + fcy*fny + fcz*fnz)
          END DO ! ifcl      
       END IF ! v4
    END IF ! global%solverType
    
! ==============================================================================
!   Normalize grid speed
! ==============================================================================    
              
    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
      vsnew = THRD*vsnew
      vsold = THRD*vsold

      pGrid%gs(ifc) = (3.0_RFREAL*vsnew - vsold)/ &
                      (2.0_RFREAL*pGrid%fn(XYZMAG,ifc)*global%dtMin)
    ELSE
      pGrid%gs(ifc) = THRD*pGrid%gs(ifc)/(pGrid%fn(XYZMAG,ifc)*global%dtMin) 
    END IF ! global%solverType
  END DO ! ifc

! ******************************************************************************
! Write out extrema for interior faces
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    IF ( pGrid%nFaces > 0 ) THEN ! Have gm test case with 0 faces...
      WRITE(STDOUT,'(A,5X,A,2(1X,E24.15),2(1X,I9))') & 
            SOLVER_NAME,'Interior:', & 
            MINVAL(pGrid%gs),MAXVAL(pGrid%gs), & 
            MINLOC(pGrid%gs),MAXLOC(pGrid%gs)
    END IF ! pGrid%nFacesTot    
  END IF ! global%myProcid  

! ******************************************************************************
! Boundary faces
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
        
! ==============================================================================
!   Moving faces
! ==============================================================================        
        
    IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN                 
      DO ifc = 1,pPatch%nBFaces
        v1 = pPatch%bv(pPatch%bf2v(1,ifc))
        v2 = pPatch%bv(pPatch%bf2v(2,ifc))
        v3 = pPatch%bv(pPatch%bf2v(3,ifc))

        pPatch%gs(ifc) = 0.0_RFREAL

        IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
           vsnew = 0.0_RFREAL
           vsold = 0.0_RFREAL
        END IF ! global%solverType
      
! ------------------------------------------------------------------------------
!     Triangular face, sweeps out prism during grid motion
! ------------------------------------------------------------------------------
            
        IF ( pPatch%bf2v(4,ifc) == VERT_NONE ) THEN ! triangular face
          xyzPri(XCOORD,1) = pXyzOld(XCOORD,v1)
          xyzPri(XCOORD,2) = pXyzOld(XCOORD,v2)
          xyzPri(XCOORD,3) = pXyzOld(XCOORD,v3)            

          xyzPri(YCOORD,1) = pXyzOld(YCOORD,v1)
          xyzPri(YCOORD,2) = pXyzOld(YCOORD,v2)
          xyzPri(YCOORD,3) = pXyzOld(YCOORD,v3)     

          xyzPri(ZCOORD,1) = pXyzOld(ZCOORD,v1)
          xyzPri(ZCOORD,2) = pXyzOld(ZCOORD,v2)
          xyzPri(ZCOORD,3) = pXyzOld(ZCOORD,v3)      

          xyzPri(XCOORD,4) = pXyz(XCOORD,v1)
          xyzPri(XCOORD,5) = pXyz(XCOORD,v2)
          xyzPri(XCOORD,6) = pXyz(XCOORD,v3)            

          xyzPri(YCOORD,4) = pXyz(YCOORD,v1)
          xyzPri(YCOORD,5) = pXyz(YCOORD,v2)
          xyzPri(YCOORD,6) = pXyz(YCOORD,v3)     

          xyzPri(ZCOORD,4) = pXyz(ZCOORD,v1)
          xyzPri(ZCOORD,5) = pXyz(ZCOORD,v2)
          xyzPri(ZCOORD,6) = pXyz(ZCOORD,v3)      

! ------- Compute average coordinates (approximate centroid) 

          xyzAvg(XCOORD) = 0.0_RFREAL
          xyzAvg(YCOORD) = 0.0_RFREAL
          xyzAvg(ZCOORD) = 0.0_RFREAL            

          term = 1.0_RFREAL/6.0_RFREAL

          DO iv = 1,6
            xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzPri(XCOORD,iv)
            xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzPri(YCOORD,iv)
            xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzPri(ZCOORD,iv)                
          END DO ! iv          
        
          xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
          xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
          xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                   
        
! ------- Loop over the five faces of the prism to compute volume
    
          DO ifcl = 1,5
            v1l = f2vPri(1,ifcl)
            v2l = f2vPri(2,ifcl)
            v3l = f2vPri(3,ifcl)
            v4l = f2vPri(4,ifcl)

            IF ( v4l == 0 ) THEN ! triangular face
              xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
              xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
              xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)

              xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
              xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
              xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)

              xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
              xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
              xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)

              CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fnx,fny,fnz)
              CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fcx,fcy,fcz)        
            ELSE 
              xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
              xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
              xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
              xyzNodes(XCOORD,4) = xyzPri(XCOORD,v4l) - xyzAvg(XCOORD)          

              xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
              xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
              xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
              xyzNodes(YCOORD,4) = xyzPri(YCOORD,v4l) - xyzAvg(YCOORD)

              xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
              xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
              xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
              xyzNodes(ZCOORD,4) = xyzPri(ZCOORD,v4l) - xyzAvg(ZCOORD)                

              CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
              CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
            END IF ! v4l

            IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
               vsnew = vsnew + (fcx*fnx + fcy*fny + fcz*fnz)
            ELSE
               pPatch%gs(ifc) = pPatch%gs(ifc) + (fcx*fnx + fcy*fny + fcz*fnz)
            END IF ! global%solverType

          END DO ! ifcl      

! ------------------------------------------------------------------------------
!     Quadrilateral face, sweeps out hexahedron during grid motion
! ------------------------------------------------------------------------------
      
        ELSE ! quadrilateral face
          v4 = pPatch%bv(pPatch%bf2v(4,ifc))

          xyzHex(XCOORD,1) = pXyzOld(XCOORD,v1)
          xyzHex(XCOORD,2) = pXyzOld(XCOORD,v2)
          xyzHex(XCOORD,3) = pXyzOld(XCOORD,v3)            
          xyzHex(XCOORD,4) = pXyzOld(XCOORD,v4)    

          xyzHex(YCOORD,1) = pXyzOld(YCOORD,v1)
          xyzHex(YCOORD,2) = pXyzOld(YCOORD,v2)
          xyzHex(YCOORD,3) = pXyzOld(YCOORD,v3)     
          xyzHex(YCOORD,4) = pXyzOld(YCOORD,v4)  

          xyzHex(ZCOORD,1) = pXyzOld(ZCOORD,v1)
          xyzHex(ZCOORD,2) = pXyzOld(ZCOORD,v2)
          xyzHex(ZCOORD,3) = pXyzOld(ZCOORD,v3)
          xyzHex(ZCOORD,4) = pXyzOld(ZCOORD,v4)              

          xyzHex(XCOORD,5) = pXyz(XCOORD,v1)
          xyzHex(XCOORD,6) = pXyz(XCOORD,v2)
          xyzHex(XCOORD,7) = pXyz(XCOORD,v3) 
          xyzHex(XCOORD,8) = pXyz(XCOORD,v4)                  

          xyzHex(YCOORD,5) = pXyz(YCOORD,v1)
          xyzHex(YCOORD,6) = pXyz(YCOORD,v2)
          xyzHex(YCOORD,7) = pXyz(YCOORD,v3)     
          xyzHex(YCOORD,8) = pXyz(YCOORD,v4)   

          xyzHex(ZCOORD,5) = pXyz(ZCOORD,v1)
          xyzHex(ZCOORD,6) = pXyz(ZCOORD,v2)
          xyzHex(ZCOORD,7) = pXyz(ZCOORD,v3)   
          xyzHex(ZCOORD,8) = pXyz(ZCOORD,v4)            

! ------- Compute average coordinates (approximate centroid)

          xyzAvg(XCOORD) = 0.0_RFREAL
          xyzAvg(YCOORD) = 0.0_RFREAL
          xyzAvg(ZCOORD) = 0.0_RFREAL            

          term = 1.0_RFREAL/8.0_RFREAL

          DO iv = 1,8
            xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzHex(XCOORD,iv)
            xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzHex(YCOORD,iv)
            xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzHex(ZCOORD,iv)                
          END DO ! iv       

          xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
          xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
          xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                   

!-------- Loop over the six faces of the hexahedron to compute volume

          DO ifcl = 1,6
            v1l = f2vHex(1,ifcl)
            v2l = f2vHex(2,ifcl)
            v3l = f2vHex(3,ifcl)
            v4l = f2vHex(4,ifcl)

            xyzNodes(XCOORD,1) = xyzHex(XCOORD,v1l) - xyzAvg(XCOORD)
            xyzNodes(XCOORD,2) = xyzHex(XCOORD,v2l) - xyzAvg(XCOORD)
            xyzNodes(XCOORD,3) = xyzHex(XCOORD,v3l) - xyzAvg(XCOORD)
            xyzNodes(XCOORD,4) = xyzHex(XCOORD,v4l) - xyzAvg(XCOORD)       

            xyzNodes(YCOORD,1) = xyzHex(YCOORD,v1l) - xyzAvg(YCOORD)
            xyzNodes(YCOORD,2) = xyzHex(YCOORD,v2l) - xyzAvg(YCOORD)
            xyzNodes(YCOORD,3) = xyzHex(YCOORD,v3l) - xyzAvg(YCOORD)
            xyzNodes(YCOORD,4) = xyzHex(YCOORD,v4l) - xyzAvg(YCOORD)       

            xyzNodes(ZCOORD,1) = xyzHex(ZCOORD,v1l) - xyzAvg(ZCOORD)
            xyzNodes(ZCOORD,2) = xyzHex(ZCOORD,v2l) - xyzAvg(ZCOORD)
            xyzNodes(ZCOORD,3) = xyzHex(ZCOORD,v3l) - xyzAvg(ZCOORD)
            xyzNodes(ZCOORD,4) = xyzHex(ZCOORD,v4l) - xyzAvg(ZCOORD)       

            CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
            CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 

            IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
               vsnew = vsnew + (fcx*fnx + fcy*fny + fcz*fnz)
            ELSE
               pPatch%gs(ifc) = pPatch%gs(ifc) + (fcx*fnx + fcy*fny + fcz*fnz)
            END IF ! global%solverType

          END DO ! ifcl            
        END IF ! v4

! ******************************************************************************
! Boundary faces with Old2 grid for implicit code
! ******************************************************************************

        IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN

! ------------------------------------------------------------------------------
!     Triangular face, sweeps out prism during grid motion
! ------------------------------------------------------------------------------
            
           IF ( pPatch%bf2v(4,ifc) == VERT_NONE ) THEN ! triangular face
              xyzPri(XCOORD,1) = pXyzOld2(XCOORD,v1)
              xyzPri(XCOORD,2) = pXyzOld2(XCOORD,v2)
              xyzPri(XCOORD,3) = pXyzOld2(XCOORD,v3)            

              xyzPri(YCOORD,1) = pXyzOld2(YCOORD,v1)
              xyzPri(YCOORD,2) = pXyzOld2(YCOORD,v2)
              xyzPri(YCOORD,3) = pXyzOld2(YCOORD,v3)     

              xyzPri(ZCOORD,1) = pXyzOld2(ZCOORD,v1)
              xyzPri(ZCOORD,2) = pXyzOld2(ZCOORD,v2)
              xyzPri(ZCOORD,3) = pXyzOld2(ZCOORD,v3)      
              
              xyzPri(XCOORD,4) = pXyzOld(XCOORD,v1)
              xyzPri(XCOORD,5) = pXyzOld(XCOORD,v2)
              xyzPri(XCOORD,6) = pXyzOld(XCOORD,v3)            
              
              xyzPri(YCOORD,4) = pXyzOld(YCOORD,v1)
              xyzPri(YCOORD,5) = pXyzOld(YCOORD,v2)
              xyzPri(YCOORD,6) = pXyzOld(YCOORD,v3)     
              
              xyzPri(ZCOORD,4) = pXyzOld(ZCOORD,v1)
              xyzPri(ZCOORD,5) = pXyzOld(ZCOORD,v2)
              xyzPri(ZCOORD,6) = pXyzOld(ZCOORD,v3)      

! ------- Compute average coordinates (approximate centroid) 

              xyzAvg(XCOORD) = 0.0_RFREAL
              xyzAvg(YCOORD) = 0.0_RFREAL
              xyzAvg(ZCOORD) = 0.0_RFREAL            
              
              term = 1.0_RFREAL/6.0_RFREAL
              
              DO iv = 1,6
                 xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzPri(XCOORD,iv)
                 xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzPri(YCOORD,iv)
                 xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzPri(ZCOORD,iv)                
              END DO ! iv          
              
              xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
              xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
              xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                   
        
! ------- Loop over the five faces of the prism to compute volume
    
              DO ifcl = 1,5
                 v1l = f2vPri(1,ifcl)
                 v2l = f2vPri(2,ifcl)
                 v3l = f2vPri(3,ifcl)
                 v4l = f2vPri(4,ifcl)
                 
                 IF ( v4l == 0 ) THEN ! triangular face
                    xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
                    xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
                    xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
                    
                    xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
                    xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
                    xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
                    
                    xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
                    xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
                    xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
                    
                    CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fnx,fny,fnz)
                    CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fcx,fcy,fcz)        
                 ELSE 
                    xyzNodes(XCOORD,1) = xyzPri(XCOORD,v1l) - xyzAvg(XCOORD)
                    xyzNodes(XCOORD,2) = xyzPri(XCOORD,v2l) - xyzAvg(XCOORD)
                    xyzNodes(XCOORD,3) = xyzPri(XCOORD,v3l) - xyzAvg(XCOORD)
                    xyzNodes(XCOORD,4) = xyzPri(XCOORD,v4l) - xyzAvg(XCOORD)          
                    
                    xyzNodes(YCOORD,1) = xyzPri(YCOORD,v1l) - xyzAvg(YCOORD)
                    xyzNodes(YCOORD,2) = xyzPri(YCOORD,v2l) - xyzAvg(YCOORD)
                    xyzNodes(YCOORD,3) = xyzPri(YCOORD,v3l) - xyzAvg(YCOORD)
                    xyzNodes(YCOORD,4) = xyzPri(YCOORD,v4l) - xyzAvg(YCOORD)
                    
                    xyzNodes(ZCOORD,1) = xyzPri(ZCOORD,v1l) - xyzAvg(ZCOORD)
                    xyzNodes(ZCOORD,2) = xyzPri(ZCOORD,v2l) - xyzAvg(ZCOORD)
                    xyzNodes(ZCOORD,3) = xyzPri(ZCOORD,v3l) - xyzAvg(ZCOORD)
                    xyzNodes(ZCOORD,4) = xyzPri(ZCOORD,v4l) - xyzAvg(ZCOORD)                
                    
                    CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
                    CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
                 END IF ! v4l
                 
                 vsold = vsold + (fcx*fnx + fcy*fny + fcz*fnz)
                 
              END DO ! ifcl      

! ------------------------------------------------------------------------------
!     Quadrilateral face, sweeps out hexahedron during grid motion
! ------------------------------------------------------------------------------
      
           ELSE ! quadrilateral face
              v4 = pPatch%bv(pPatch%bf2v(4,ifc))
              
              xyzHex(XCOORD,1) = pXyzOld(XCOORD,v1)
              xyzHex(XCOORD,2) = pXyzOld(XCOORD,v2)
              xyzHex(XCOORD,3) = pXyzOld(XCOORD,v3)            
              xyzHex(XCOORD,4) = pXyzOld(XCOORD,v4)    
              
              xyzHex(YCOORD,1) = pXyzOld(YCOORD,v1)
              xyzHex(YCOORD,2) = pXyzOld(YCOORD,v2)
              xyzHex(YCOORD,3) = pXyzOld(YCOORD,v3)     
              xyzHex(YCOORD,4) = pXyzOld(YCOORD,v4)  
              
              xyzHex(ZCOORD,1) = pXyzOld(ZCOORD,v1)
              xyzHex(ZCOORD,2) = pXyzOld(ZCOORD,v2)
              xyzHex(ZCOORD,3) = pXyzOld(ZCOORD,v3)
              xyzHex(ZCOORD,4) = pXyzOld(ZCOORD,v4)              
              
              xyzHex(XCOORD,5) = pXyz(XCOORD,v1)
              xyzHex(XCOORD,6) = pXyz(XCOORD,v2)
              xyzHex(XCOORD,7) = pXyz(XCOORD,v3) 
              xyzHex(XCOORD,8) = pXyz(XCOORD,v4)                  
              
              xyzHex(YCOORD,5) = pXyz(YCOORD,v1)
              xyzHex(YCOORD,6) = pXyz(YCOORD,v2)
              xyzHex(YCOORD,7) = pXyz(YCOORD,v3)     
              xyzHex(YCOORD,8) = pXyz(YCOORD,v4)   
              
              xyzHex(ZCOORD,5) = pXyz(ZCOORD,v1)
              xyzHex(ZCOORD,6) = pXyz(ZCOORD,v2)
              xyzHex(ZCOORD,7) = pXyz(ZCOORD,v3)   
              xyzHex(ZCOORD,8) = pXyz(ZCOORD,v4)            

! ------- Compute average coordinates (approximate centroid)

              xyzAvg(XCOORD) = 0.0_RFREAL
              xyzAvg(YCOORD) = 0.0_RFREAL
              xyzAvg(ZCOORD) = 0.0_RFREAL            
              
              term = 1.0_RFREAL/8.0_RFREAL
              
              DO iv = 1,8
                 xyzAvg(XCOORD) = xyzAvg(XCOORD) + xyzHex(XCOORD,iv)
                 xyzAvg(YCOORD) = xyzAvg(YCOORD) + xyzHex(YCOORD,iv)
                 xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + xyzHex(ZCOORD,iv)                
              END DO ! iv       
              
              xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
              xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
              xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)                   
              
!-------- Loop over the six faces of the hexahedron to compute volume

              DO ifcl = 1,6
                 v1l = f2vHex(1,ifcl)
                 v2l = f2vHex(2,ifcl)
                 v3l = f2vHex(3,ifcl)
                 v4l = f2vHex(4,ifcl)
                 
                 xyzNodes(XCOORD,1) = xyzHex(XCOORD,v1l) - xyzAvg(XCOORD)
                 xyzNodes(XCOORD,2) = xyzHex(XCOORD,v2l) - xyzAvg(XCOORD)
                 xyzNodes(XCOORD,3) = xyzHex(XCOORD,v3l) - xyzAvg(XCOORD)
                 xyzNodes(XCOORD,4) = xyzHex(XCOORD,v4l) - xyzAvg(XCOORD)       
                 
                 xyzNodes(YCOORD,1) = xyzHex(YCOORD,v1l) - xyzAvg(YCOORD)
                 xyzNodes(YCOORD,2) = xyzHex(YCOORD,v2l) - xyzAvg(YCOORD)
                 xyzNodes(YCOORD,3) = xyzHex(YCOORD,v3l) - xyzAvg(YCOORD)
                 xyzNodes(YCOORD,4) = xyzHex(YCOORD,v4l) - xyzAvg(YCOORD)       
                 
                 xyzNodes(ZCOORD,1) = xyzHex(ZCOORD,v1l) - xyzAvg(ZCOORD)
                 xyzNodes(ZCOORD,2) = xyzHex(ZCOORD,v2l) - xyzAvg(ZCOORD)
                 xyzNodes(ZCOORD,3) = xyzHex(ZCOORD,v3l) - xyzAvg(ZCOORD)
                 xyzNodes(ZCOORD,4) = xyzHex(ZCOORD,v4l) - xyzAvg(ZCOORD)       
                 
                 CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fnx,fny,fnz) 
                 CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fcx,fcy,fcz) 
                 
                 vsold = vsold + (fcx*fnx + fcy*fny + fcz*fnz)
                 
              END DO ! ifcl            
           END IF ! v4
           
        END IF ! global%solverType

! ------------------------------------------------------------------------------
!       Normalize grid speed
! ------------------------------------------------------------------------------

        IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
           vsnew = THRD*vsnew
           vsold = THRD*vsold
           pPatch%gs(ifc) = (3.0_RFREAL*vsnew - vsold)/(2.0_RFREAL*pPatch%fn(XYZMAG,ifc)*global%dtMin)
        ELSE
           pPatch%gs(ifc) = THRD*pPatch%gs(ifc)/ (pPatch%fn(XYZMAG,ifc)*global%dtMin)
        END IF ! global%solverType
    
    END DO ! ifc
 
! ==============================================================================
!   Non-moving patches (should have precisely zero grid speed, so enforce)
! ============================================================================== 
      
    ELSE 
      DO ifc = 1,pPatch%nBFaces
        pPatch%gs(ifc) = 0.0_RFREAL
      END DO ! ifc
    END IF ! pPatch%movePatchDir    
               
! ==============================================================================
!   Write out extrema for patches with actual faces
! ==============================================================================     
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      IF ( pPatch%nBFaces > 0 ) THEN ! can have zero faces   
        WRITE(STDOUT,'(A,5X,A,I3,A,2(1X,E24.15),2(1X,I9))') & 
              SOLVER_NAME,'Patch',iPatch,':', & 
              MINVAL(pPatch%gs(1:pPatch%nBFaces)), &
              MAXVAL(pPatch%gs(1:pPatch%nBFaces)), & 
              MINLOC(pPatch%gs(1:pPatch%nBFaces)), &
              MAXLOC(pPatch%gs(1:pPatch%nBFaces)) 
      END IF ! pPatch%nBFaces
    END IF ! global%myProcid         
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid speeds done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)


END SUBROUTINE RFLU_ComputeGridSpeeds

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeGridSpeeds.F90,v $
! Revision 1.9  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/02/08 21:31:48  hdewey2
! Added grid speed computation for implicit solver
!
! Revision 1.6  2005/06/09 20:22:58  haselbac
! Replaced movePatch by movePatchDir
!
! Revision 1.5  2004/10/19 19:41:18  haselbac
! Cosmetics only
!
! Revision 1.4  2003/03/15 18:28:01  haselbac
! Changed loop limits, added VERT_NONE, other changes for || gm
!
! Revision 1.3  2003/01/28 14:31:27  haselbac
! Set gs to zero on non-moving boundaries, use scaled coordinates, clean-up
!
! Revision 1.2  2002/11/08 21:29:45  haselbac
! Some cosmetics and clean-up
!
! Revision 1.1  2002/10/27 19:20:08  haselbac
! Initial revision
!
! ******************************************************************************







