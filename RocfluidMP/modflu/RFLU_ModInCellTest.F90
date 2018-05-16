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
! Purpose: Suite of routines related to in-cell test.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModInCellTest.F90,v 1.4 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModInCellTest

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid    
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_ICT_ComputeTolerance, & 
            RFLU_ICT_TestFacePlanar, & 
            RFLU_ICT_TestFaceQuadBilinear, & 
            RFLU_ICT_TestInCell, & 
            RFLU_ICT_TestInCellFancy, & 
            RFLU_ICT_TestInCellLohner
      
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModInCellTest.F90,v $ $Revision: 1.4 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! *****************************************************************************
!
! Purpose: Compute tolerance for quadrilateral faces for in-cell test. 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None. 
!
! *****************************************************************************

SUBROUTINE RFLU_ICT_ComputeTolerance(pRegion)
   
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: ifg,ifl,iPatch,ivg,ivl
  REAL(RFREAL) :: eps,epsMax,nx,ny,nz,safetyFactor,xc,yc,zc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_ComputeTolerance',&
  'RFLU_ModInCellTest.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing tolerance for in-cell test...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal
  END IF ! global%myProcid

! =============================================================================
! Set grid pointer and initialize variables
! =============================================================================

  pGrid => pRegion%grid

  epsMax = -HUGE(1.0_RFREAL)

  safetyFactor = 2.0_RFREAL

! *****************************************************************************
! Compute plane equation defect 
! *****************************************************************************

  DO ifg = 1,pGrid%nFaces
    IF ( pGrid%f2v(4,ifg) /= VERT_NONE ) THEN 
      xc = pGrid%fc(XCOORD,ifg)
      yc = pGrid%fc(YCOORD,ifg)
      zc = pGrid%fc(ZCOORD,ifg)
      
      nx = pGrid%fn(XCOORD,ifg)
      ny = pGrid%fn(YCOORD,ifg)
      nz = pGrid%fn(ZCOORD,ifg)

      DO ivl = 1,4
        ivg = pGrid%f2v(ivl,ifg)

        eps = (pGrid%xyz(XCOORD,ivg) - xc)*nx & 
            + (pGrid%xyz(YCOORD,ivg) - yc)*ny & 
            + (pGrid%xyz(ZCOORD,ivg) - zc)*nz

        epsMax = MAX(eps,epsMax)
      END DO ! ivl
    END IF ! pGrid%f2v
  END DO ! ifg

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    DO ifl = 1,pPatch%nBFaces
      IF ( pPatch%bf2v(4,ifl) /= VERT_NONE ) THEN 
        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)
      
        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)

        DO ivl = 1,4
          ivg = pPatch%bv(pPatch%bf2v(ivl,ifl))
  
          eps = (pGrid%xyz(XCOORD,ivg) - xc)*nx & 
              + (pGrid%xyz(YCOORD,ivg) - yc)*ny & 
              + (pGrid%xyz(ZCOORD,ivg) - zc)*nz

          epsMax = MAX(eps,epsMax)
        END DO ! ivl
      END IF ! pPatch%bf2v
    END DO ! ifl
  END DO ! iPatch

! *****************************************************************************
! Write info
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,E13.6)') SOLVER_NAME, & 
                                      'Maximum face planarity defect:',epsMax 
  END IF ! global%myProcid

! *****************************************************************************
! Set new tolerance
! *****************************************************************************

  IF ( epsMax > pRegion%mixtInput%tolerICT ) THEN 
    pRegion%mixtInput%tolerICT = safetyFactor*epsMax

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,12X,E13.6)') SOLVER_NAME,'Tolerance reset to:', &
                                         pRegion%mixtInput%tolerICT
    END IF ! global%myProcid
  END IF ! epsMax 

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Computing tolerance for in-cell test done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ICT_ComputeTolerance







! *****************************************************************************
!
! Purpose: Compute contribution to in-cell test from planar triangular or
!   quadrilateral face. 
!
! Description: Compute dot product of face normal vector and relative vector
!   between face centroid and location in question.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   icg                 Global cell index
!   iPatch              Patch index
!   ifg                 Face index
!
! Output: 
!   dotp                Dot product
!
! Notes: 
!   1. Only applicable to planar faces.
!
! *****************************************************************************

SUBROUTINE RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)
   
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icg,ifg,iPatch
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  REAL(RFREAL), INTENT(OUT) :: dotp
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: c1,c2
  REAL(RFREAL) :: fnx,fny,fnz,xCofg,yCofg,zCofg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_TestFacePlanar',&
  'RFLU_ModInCellTest.F90')

! =============================================================================
! Set grid pointer
! =============================================================================

  pGrid => pRegion%grid

! *****************************************************************************
! Compute dot product
! *****************************************************************************

  IF ( iPatch == 0 ) THEN ! interior face
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    xCofg = pGrid%fc(XCOORD,ifg)
    yCofg = pGrid%fc(YCOORD,ifg)
    zCofg = pGrid%fc(ZCOORD,ifg)                    

    IF ( c1 == icg ) THEN 
      fnx =  pGrid%fn(XCOORD,ifg)
      fny =  pGrid%fn(YCOORD,ifg)
      fnz =  pGrid%fn(ZCOORD,ifg)
    ELSE IF ( c2 == icg ) THEN
      fnx = -pGrid%fn(XCOORD,ifg)
      fny = -pGrid%fn(YCOORD,ifg)
      fnz = -pGrid%fn(ZCOORD,ifg)
    ELSE ! defensive programming
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! c1                    
  ELSE ! boundary face  
    pPatch => pRegion%patches(iPatch)

    xCofg = pPatch%fc(XCOORD,ifg)
    yCofg = pPatch%fc(YCOORD,ifg)
    zCofg = pPatch%fc(ZCOORD,ifg)                    

    fnx = pPatch%fn(XCOORD,ifg)
    fny = pPatch%fn(YCOORD,ifg)
    fnz = pPatch%fn(ZCOORD,ifg)          
  END IF ! iPatch

  dotp = (xCofg-xLoc)*fnx + (yCofg-yLoc)*fny + (zCofg-zLoc)*fnz

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ICT_TestFacePlanar






! ******************************************************************************
!
! Purpose: Compute contribution to in-cell test from bilinear quadrilateral 
!   face. 
!
! Description: A point is located inside a cell if the dot products of the 
!   outward-facing face normal vectors and the relative vector of the face-
!   centroids and the point in question are all positive.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   iPatch              Patch index
!   ifg                 Face index
!
! Output: 
!   dotp                Dot product indicating whether point inside cell
!
! Notes: 
!   1. Return from this routine early if detect negative dot product.
!
! ******************************************************************************

SUBROUTINE RFLU_ICT_TestFaceQuadBilinear(pRegion,xLoc,yLoc,zLoc,iPatch,ifg,dotp)
      
  USE RFLU_ModBilinearPatch    
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ifg,iPatch
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  REAL(RFREAL), INTENT(OUT) :: dotp
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: v1,v2,v3,v4
  REAL(RFREAL) :: ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,nx,ny,nz,u,v,x,x1,x2, & 
                  x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_TestFaceQuadBilinear',&
  'RFLU_ModInCellTest.F90')

! ******************************************************************************
! Set pointer
! ******************************************************************************

  pGrid => pRegion%grid
  
! ******************************************************************************
! Get vertices
! ******************************************************************************

  IF ( iPatch == 0 ) THEN ! interior face
    v1 = pGrid%f2v(1,ifg)
    v2 = pGrid%f2v(2,ifg)
    v3 = pGrid%f2v(3,ifg)
    v4 = pGrid%f2v(4,ifg)    

    IF ( v4 == VERT_NONE ) THEN ! defensive coding
! TEMPORARY    
      WRITE(*,*) 'ERROR!'
      STOP
! END TEMPORARY
    END IF ! v4                                   
  ELSE ! boundary face  
    pPatch => pRegion%patches(iPatch)

    v1 = pPatch%bv(pPatch%bf2v(1,ifg))
    v2 = pPatch%bv(pPatch%bf2v(2,ifg))
    v3 = pPatch%bv(pPatch%bf2v(3,ifg))
                     
    IF ( pPatch%bf2v(4,ifg) == VERT_NONE ) THEN ! defensive coding
! TEMPORARY    
      WRITE(*,*) 'ERROR!'
      STOP
! END TEMPORARY
    END IF ! pPatch%bf2v(4,ifg)
    
    v4 = pPatch%bv(pPatch%bf2v(4,ifg))              
  END IF ! iloc

  x1 = pGrid%xyz(XCOORD,v1)
  x2 = pGrid%xyz(XCOORD,v2)
  x3 = pGrid%xyz(XCOORD,v3)   
  x4 = pGrid%xyz(XCOORD,v4) 

  y1 = pGrid%xyz(YCOORD,v1)
  y2 = pGrid%xyz(YCOORD,v2)
  y3 = pGrid%xyz(YCOORD,v3)   
  y4 = pGrid%xyz(YCOORD,v4) 

  z1 = pGrid%xyz(ZCOORD,v1)
  z2 = pGrid%xyz(ZCOORD,v2)
  z3 = pGrid%xyz(ZCOORD,v3)   
  z4 = pGrid%xyz(ZCOORD,v4) 

! ******************************************************************************
! Compute geometric terms of parametric patch representation
! ******************************************************************************

  ax = x1 - x2 - x4 + x3
  ay = y1 - y2 - y4 + y3
  az = z1 - z2 - z4 + z3 

  bx = x2 - x1
  by = y2 - y1
  bz = z2 - z1  

  cx = x4 - x1
  cy = y4 - y1
  cz = z4 - z1  

  dx = x1
  dy = y1
  dz = z1
  
! ******************************************************************************
! Compute closest point on bilinear patch 
! ******************************************************************************

  CALL RFLU_BLIN_FindClosestPoint(global,ax,ay,az,bx,by,bz,cx,cy,cz, & 
                                  dx,dy,dz,xLoc,yLoc,zLoc,u,v,x,y,z)  
   
! ******************************************************************************
! Compute intersection distance (here referred to simply as dot product)
! ******************************************************************************

  IF ( (u < 0.0_RFREAL .OR. u > 1.0_RFREAL) .OR. & 
       (v < 0.0_RFREAL .OR. v > 1.0_RFREAL) ) THEN 
    dotp = CRAZY_VALUE_INT ! Must be large negative number
  ELSE 
    CALL RFLU_BLIN_ComputeNormal(global,ax,ay,az,bx,by,bz,cx,cy,cz, &
                                 dx,dy,dz,u,v,nx,ny,nz)  

    dotp = (x-xLoc)*nx + (y-yLoc)*ny + (z-zLoc)*nz
  END IF ! u

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ICT_TestFaceQuadBilinear






! ******************************************************************************
!
! Purpose: Determine whether point is inside cell.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   icg                 Global cell index
!
! Output: 
!   RFLU_ICT_TestInCell = .TRUE.     if point inside cell
!   RFLU_ICT_TestInCell = .FALSE.    if point outside cell
!
! Notes: 
!   1. Return from this routine early if detect negative dot product.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg)
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: cntr,icl,ict,ifg,ifl,iPatch,nFaces,v4
  REAL(RFREAL) :: dotp,toler
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_TestInCell',&
  'RFLU_ModInCellTest.F90')

! ==============================================================================
! Set grid pointer and initialize
! ==============================================================================

  pGrid => pRegion%grid
  
  RFLU_ICT_TestInCell = .FALSE.

  toler = -pRegion%mixtInput%tolerICT 

! ******************************************************************************
! Check whether point is inside cell
! ******************************************************************************

  ict = pGrid%cellGlob2Loc(1,icg) ! cell type
  icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

  cntr = 0

  SELECT CASE ( ict )

! ==============================================================================
!   Tetrahedron
! ==============================================================================

    CASE ( CELL_TYPE_TET )
      nFaces = 4
    
      DO ifl = 1,nFaces
        iPatch = pGrid%tet2f(1,ifl,icl)
        ifg    = pGrid%tet2f(2,ifl,icl)        

        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)

        IF ( dotp >= toler ) THEN 
          cntr = cntr + 1
        ELSE 
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifc     

! ==============================================================================
!   Hexahedron
! ==============================================================================

    CASE ( CELL_TYPE_HEX ) 
      nFaces = 6
    
      DO ifl = 1,nFaces
        iPatch = pGrid%hex2f(1,ifl,icl)
        ifg    = pGrid%hex2f(2,ifl,icl)        

! TEMPORARY
        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)
!        CALL RFLU_ICT_TestFaceQuadBilinear(pRegion,xLoc,yLoc,zLoc,iPatch,ifg, &
!                                           dotp)
! END TEMPORARY

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifc   

! ==============================================================================
!   Prism
! ==============================================================================

    CASE ( CELL_TYPE_PRI )                 
      nFaces = 5
    
      DO ifl = 1,nFaces
        iPatch = pGrid%pri2f(1,ifl,icl)
        ifg    = pGrid%pri2f(2,ifl,icl)        

        IF ( iPatch == 0 ) THEN 
          v4 = pGrid%f2v(4,ifg)
        ELSE IF ( iPatch > 0 ) THEN
          pPatch => pRegion%patches(iPatch)
         
          v4 = pPatch%bf2v(4,ifg)
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! iPatch        

        IF ( v4 == VERT_NONE ) THEN 
          CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg, &
                                       dotp)  
        ELSE
! TEMPORARY
          CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg, &
                                       dotp)  
!          CALL RFLU_ICT_TestFaceQuadBilinear(pRegion,xLoc,yLoc,zLoc,iPatch, &
!                                             ifg,dotp)
! END TEMPORARY
        END IF ! v4

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifc 

! ==============================================================================
!   Pyramid
! ==============================================================================

    CASE ( CELL_TYPE_PYR ) 
      nFaces = 5
    
      DO ifl = 1,nFaces
        iPatch = pGrid%pyr2f(1,ifl,icl)
        ifg    = pGrid%pyr2f(2,ifl,icl)        

        IF ( iPatch == 0 ) THEN 
          v4 = pGrid%f2v(4,ifg)
        ELSE IF ( iPatch > 0 ) THEN
          pPatch => pRegion%patches(iPatch)
         
          v4 = pPatch%bf2v(4,ifg)
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! iPatch        

        IF ( v4 == VERT_NONE ) THEN 
          CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg, &
                                       dotp)  
        ELSE
! TEMPORARY
          CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg, &
                                       dotp)  
!          CALL RFLU_ICT_TestFaceQuadBilinear(pRegion,xLoc,yLoc,zLoc,iPatch, &
!                                             ifg,dotp)
! END TEMPORARY
        END IF ! v4
        
        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifc 

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                      
  END SELECT ! ict        

! ==============================================================================
! Set RFLU_ICT_TestInCell to TRUE if have nFaces positive dot products
! ==============================================================================

  IF ( cntr == nFaces ) THEN 
    RFLU_ICT_TestInCell = .TRUE. 
  END IF ! cntr

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_ICT_TestInCell







! *****************************************************************************
!
! Purpose: Determine whether point is inside cell and return information about
!   which face which might have lead to failure.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   icg                 Global cell index
!
! Output: 
!   testInCell          True if location in cell, false otherwise
!   iPatchOut           Patch index of face which failed test
!   ifgOut              Global index of face which failed test
!
! Notes: 
!   1. Return from this routine early if detect negative dot product.
!
! *****************************************************************************

SUBROUTINE RFLU_ICT_TestInCellFancy(pRegion,xLoc,yLoc,zLoc,icg,testInCell, & 
                                    iPatchOut,ifgOut)
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  LOGICAL, INTENT(OUT) :: testInCell
  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(OUT) :: iPatchOut,ifgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: cntr,icl,ict,ifg,ifl,iPatch,nFaces
  REAL(RFREAL) :: dotp,toler
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_TestInCellFancy',&
  'RFLU_ModInCellTest.F90')

! =============================================================================
! Set grid pointer and initialize
! =============================================================================

  pGrid => pRegion%grid
  
  testInCell = .FALSE.

  toler = -pRegion%mixtInput%tolerICT

! *****************************************************************************
! Check whether point is inside cell
! *****************************************************************************

  ict = pGrid%cellGlob2Loc(1,icg) ! cell type
  icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

  cntr = 0

  SELECT CASE ( ict )
    CASE ( CELL_TYPE_TET )
      nFaces = 4
    
      DO ifl = 1,nFaces
        iPatch = pGrid%tet2f(1,ifl,icl)
        ifg    = pGrid%tet2f(2,ifl,icl)        

        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          iPatchOut = iPatch
          ifgOut    = ifg
        
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifl    
    CASE ( CELL_TYPE_HEX ) 
      nFaces = 6
    
      DO ifl = 1,nFaces
        iPatch = pGrid%hex2f(1,ifl,icl)
        ifg    = pGrid%hex2f(2,ifl,icl)        

        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          iPatchOut = iPatch
          ifgOut    = ifg

          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifl   
    CASE ( CELL_TYPE_PRI )                 
      nFaces = 5
    
      DO ifl = 1,nFaces
        iPatch = pGrid%pri2f(1,ifl,icl)
        ifg    = pGrid%pri2f(2,ifl,icl)        

        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          iPatchOut = iPatch
          ifgOut    = ifg
                
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifl 
    CASE ( CELL_TYPE_PYR ) 
      nFaces = 5
    
      DO ifl = 1,nFaces
        iPatch = pGrid%pyr2f(1,ifl,icl)
        ifg    = pGrid%pyr2f(2,ifl,icl)        

        CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg,dotp)

        IF ( dotp > toler ) THEN 
          cntr = cntr + 1
        ELSE 
          iPatchOut = iPatch
          ifgOut    = ifg
                
          CALL DeregisterFunction(global)
          RETURN
        END IF ! dotp
      END DO ! ifl 
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                      
  END SELECT          

! =============================================================================
! Set testInCell to TRUE if have nFaces positive dot products
! =============================================================================

  IF ( cntr == nFaces ) THEN 
    testInCell = .TRUE. 
    
    iPatchOut = CRAZY_VALUE_INT
    ifgOut    = CRAZY_VALUE_INT
  END IF ! cntr

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ICT_TestInCellFancy







! ******************************************************************************
!
! Purpose: Determine whether point is inside cell and return information about
!   which face failed the most, so to speak.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   icg                 Global cell index
!
! Output: 
!   testInCell          True if location in cell, false otherwise
!   ilocOut             Location of face which failed test
!   ifgOut              Global index of face which failed test
!
! Notes: 
!   1. DO NOT return from this routine early if detect negative dot product.
!   2. See R. Lohner and J. Ambrosiano, A Vectorized Particle Tracer for 
!      Unstructured Grids, J. Comp. Phys., Vol. 91, pp. 22-31, 1990.
!
! ******************************************************************************

SUBROUTINE RFLU_ICT_TestInCellLohner(pRegion,xLoc,yLoc,zLoc,icg,testInCell, & 
                                     iPatchOut,ifgOut)
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(OUT) :: testInCell
  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(OUT) :: iPatchOut,ifgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: cntr,icl,ict,ifg,ifl,iPatch,nFaces
  INTEGER, DIMENSION(1) :: iDotpMin
  INTEGER, DIMENSION(:,:), POINTER :: pC2f
  REAL(RFREAL) :: toler
  REAL(RFREAL), DIMENSION(6) :: dotp
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ICT_TestInCellLohner',&
  'RFLU_ModInCellTest.F90')

! ==============================================================================
! Set pointer and initialize 
! ==============================================================================

  pGrid => pRegion%grid
  
  testInCell = .FALSE.

  toler = -pRegion%mixtInput%tolerICT

! ******************************************************************************
! Select cell type and set pointer to cell-to-face connectivity array
! ******************************************************************************

  ict = pGrid%cellGlob2Loc(1,icg) ! cell type
  icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

  SELECT CASE ( ict ) 
    CASE ( CELL_TYPE_TET )
      pC2f => pGrid%tet2f(:,:,icl)
    CASE ( CELL_TYPE_HEX ) 
      pC2f => pGrid%hex2f(:,:,icl)
    CASE ( CELL_TYPE_PRI ) 
      pC2f => pGrid%pri2f(:,:,icl)      
    CASE ( CELL_TYPE_PYR ) 
      pC2f => pGrid%pyr2f(:,:,icl)          
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
  END SELECT ! ict
  
  nFaces = SIZE(pC2f,2)

! ******************************************************************************
! Check whether point is inside cell
! ******************************************************************************

  cntr = 0

! ==============================================================================
! Loop over faces of cell. NOTE need to initialize dotp for non-existent faces.
! ==============================================================================
    
  DO ifl = 1,nFaces
    iPatch = pC2f(1,ifl)
    ifg  = pC2f(2,ifl)        

    CALL RFLU_ICT_TestFacePlanar(pRegion,xLoc,yLoc,zLoc,icg,iPatch,ifg, &
                                 dotp(ifl))

    IF ( dotp(ifl) > toler ) THEN 
      cntr = cntr + 1
    END IF ! dotp
  END DO ! ifl   

  IF ( nFaces < 6 ) THEN 
    DO ifl = nFaces+1,6
      dotp(ifl) = HUGE(1.0_RFREAL) 
    END DO ! ifl
  END IF ! nFaces

! ==============================================================================
! Set testInCell to TRUE if have nFaces positive dot products
! ==============================================================================

  IF ( cntr == nFaces ) THEN 
    testInCell = .TRUE. 
    
    iPatchOut = CRAZY_VALUE_INT
    ifgOut    = CRAZY_VALUE_INT
  ELSE 
    iDotpMin = MINLOC(dotp(:))
  
    iPatchOut = pC2f(1,iDotpMin(1))
    ifgOut    = pC2f(2,iDotpMin(1))  
  END IF ! cntr

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ICT_TestInCellLohner





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModInCellTest


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInCellTest.F90,v $
! Revision 1.4  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.1  2005/12/24 21:17:50  haselbac
! Initial revision
!
! ******************************************************************************
  












