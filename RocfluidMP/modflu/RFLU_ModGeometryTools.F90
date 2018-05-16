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
! Purpose: Suite of geometry tools.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGeometryTools.F90,v 1.5 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2005-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGeometryTools

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_ComputeLineCellXSectFast, & 
            RFLU_ComputeLineCellXSectSafe, & 
            RFLU_TestInBoundBox, &
            RFLU_TestVectorCartAxisAligned
      
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGeometryTools.F90,v $ $Revision: 1.5 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! ******************************************************************************
!
! Purpose: Compute intersection of given line vector and faces of given cell
!   and distance between given location and intersection using fast algorithm. 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   ex,ey,ez            x-, y-, and z-components of unit line vector
!   icg                 Global cell index
!
! Output:
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of intersection
!   distOut             Distance from location to intersection 
!   iPatchOut           Face location 
!   ifgOut              Face index 
!
! Notes: 
!   1. The line vector MUST be a unit line vector. If that is not correct, the
!      distance between the given point and the intersection of the line with 
!      the faces of the given cell will not be computed correctly.
!   2. Distance might be zero if vertices lie on path, so need to include this
!      case in IF statement on dist after calls to routine which computes the
!      distance.
!   3. This is called a fast algorithm because it CANNOT detect and correct 
!      inconsistent input data, which makes it faster than the safe version.
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeLineCellXSectFast(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
                                         distOut,iPatchOut,ifgOut)
     
  USE RFLU_ModBilinearPatch, ONLY: RFLU_BLIN_ComputeXSectLine  
     
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(OUT) :: ifgOut,iPatchOut
  REAL(RFREAL), INTENT(IN) :: ex,ey,ez
  REAL(RFREAL), INTENT(INOUT) :: xLoc,yLoc,zLoc    
  REAL(RFREAL), INTENT(OUT) :: distOut
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,icl,ict,ifl,iflOut,ifg,iPatch,nFaces,nXSect
  INTEGER, DIMENSION(:,:), POINTER :: pC2f
  REAL(RFREAL) :: denom,dist,fnx,fny,fnz,numer,toler,xCofg,yCofg,zCofg
  REAL(RFREAL), DIMENSION(2) :: xsd  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLineCellXSectFast',&
  'RFLU_ModGeometryTools.F90')

! ==============================================================================
! Set grid pointer and initialize variables
! ==============================================================================

  pGrid => pRegion%grid

  toler = -pRegion%mixtInput%tolerICT ! Must be consistent with ICT tolerance

  distOut = HUGE(1.0_RFREAL)
  iflOut  = CRAZY_VALUE_INT

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
! Loop over faces of cell
! ******************************************************************************
  
  DO ifl = 1,nFaces
    iPatch = pC2f(1,ifl)
    ifg    = pC2f(2,ifl)

! ==============================================================================
!   Interior face
! ==============================================================================

    IF ( iPatch == 0 ) THEN 

! ------------------------------------------------------------------------------
!     Triangular face
! ------------------------------------------------------------------------------    
    
! TEMPORARY
!      IF ( pGrid%f2v(4,ifg) == VERT_NONE ) THEN 
! END TEMPORARY
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

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

        denom = ex*fnx + ey*fny + ez*fnz

        IF ( denom > 0.0_RFREAL ) THEN 
          xCofg = pGrid%fc(XCOORD,ifg)
          yCofg = pGrid%fc(YCOORD,ifg)
          zCofg = pGrid%fc(ZCOORD,ifg) 

          numer = (xCofg-xLoc)*fnx + (yCofg-yLoc)*fny + (zCofg-zLoc)*fnz      
          numer = MAX(numer,0.0_RFREAL) 

          dist = numer/denom

          IF ( dist < distOut ) THEN 
            distOut = dist
            iflOut  = ifl
          END IF ! dist          
        END IF ! denom 

! ------------------------------------------------------------------------------
!     Quadrilateral face
! ------------------------------------------------------------------------------    
    
! TEMPORARY
!      ELSE 
!        CALL RFLU_BLIN_ComputeXSectLine(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
!                                        iPatch,ifg,nXSect,xsd)        
!                              
!        IF ( nXSect > 0 ) THEN 
!          dist = xsd(1)
!          
!          IF ( dist < distOut ) THEN 
!            distOut = dist
!            iflOut  = ifl
!          END IF ! dist          
!        END IF ! nXSect                                   
!      END IF ! pGrid%f2v     
! END TEMPORARY

! ==============================================================================
!   Boundary face
! ==============================================================================
         
    ELSE IF ( iPatch > 0 ) THEN  
      pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!     Triangular face
! ------------------------------------------------------------------------------    

! TEMPORARY
!      IF ( pPatch%bf2v(4,ifg) == VERT_NONE ) THEN 
! END TEMPORARY
        fnx = pPatch%fn(XCOORD,ifg)
        fny = pPatch%fn(YCOORD,ifg)
        fnz = pPatch%fn(ZCOORD,ifg)          

        denom = ex*fnx + ey*fny + ez*fnz

        IF ( denom > 0.0_RFREAL ) THEN 
          xCofg = pPatch%fc(XCOORD,ifg)
          yCofg = pPatch%fc(YCOORD,ifg)
          zCofg = pPatch%fc(ZCOORD,ifg)

          numer = (xCofg-xLoc)*fnx + (yCofg-yLoc)*fny + (zCofg-zLoc)*fnz      
          numer = MAX(numer,0.0_RFREAL) 

          dist = numer/denom

          IF ( dist < distOut ) THEN 
            distOut = dist
            iflOut  = ifl
          END IF ! dist 
        END IF ! denom 

! ------------------------------------------------------------------------------
!     Quadrilateral face
! ------------------------------------------------------------------------------    

! TEMPORARY
!      ELSE 
!        CALL RFLU_BLIN_ComputeXSectLine(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
!                                        iPatch,ifg,nXSect,xsd)        
!                              
!        IF ( nXSect > 0 ) THEN 
!          dist = xsd(1)
!          
!          IF ( dist < distOut ) THEN 
!            distOut = dist
!            iflOut  = ifl
!          END IF ! dist           
!        END IF ! nXSect         
!      END IF ! pPatch%bf2v        
! END TEMPORARY

! ==============================================================================
!   Default
! ==============================================================================

    ELSE ! Defensive programming
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iPatch
  END DO ! ifl

! ******************************************************************************
! Set output
! ******************************************************************************

  iPatchOut = pC2f(1,iflOut)
  ifgOut    = pC2f(2,iflOut)      

  xLoc = xLoc + distOut*ex
  yLoc = yLoc + distOut*ey
  zLoc = zLoc + distOut*ez    

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeLineCellXSectFast








! ******************************************************************************
!
! Purpose: Compute intersection of given line vector and faces of given cell
!   and distance between given location and intersection using safe algorithm. 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   ex,ey,ez            x-, y-, and z-components of unit line vector
!   icg                 Global cell index
!
! Output:
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of intersection
!   distOut             Distance from location to intersection 
!   iPatchOut           Face location 
!   ifgOut              Face index 
!
! Notes: 
!   1. The line vector MUST be a unit line vector. If that is not correct, the
!      distance between the given point and the intersection of the line with 
!      the faces of the given cell will not be computed correctly.
!   2. Distance might be zero if vertices lie on path, so need to include this
!      case in IF statement on dist after calls to routine which computes the
!      distance.
!   3. This is called a safe algorithm because it can detect and correct 
!      inconsistent input data.    
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeLineCellXSectSafe(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
                                         distOut,iPatchOut,ifgOut)

  USE RFLU_ModBilinearPatch, ONLY: RFLU_BLIN_ComputeXSectLine   
     
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg
  INTEGER, INTENT(OUT) :: ifgOut,iPatchOut
  REAL(RFREAL), INTENT(IN) :: ex,ey,ez
  REAL(RFREAL), INTENT(INOUT) :: xLoc,yLoc,zLoc    
  REAL(RFREAL), INTENT(OUT) :: distOut
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,errorFlag,icl,ict,ifl,iflOut,iflOut1,iflOut2,ifg,iPatch, & 
             nFaces,nXSect
  INTEGER, DIMENSION(:,:), POINTER :: pC2f
  REAL(RFREAL) :: distOut1,distOut2,fnx,fny,fnz,toler,xCofg,yCofg,zCofg
  REAL(RFREAL), DIMENSION(2) :: xsd
  REAL(RFREAL), DIMENSION(6) :: denom,numer
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeLineCellXSectSafe',&
  'RFLU_ModGeometryTools.F90')

! ==============================================================================
! Set grid pointer and initialize variables
! ==============================================================================

  pGrid => pRegion%grid

  toler = -pRegion%mixtInput%tolerICT ! Must be consistent with ICT tolerance

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

! DEBUG
!  WRITE(0,*) '@@@100',pRegion%iRegionGlobal,'xsect with cell:',icg,ict,icl
! END DEBUG
                              
! ******************************************************************************
! Loop over faces of cell
! ******************************************************************************
  
  DO ifl = 1,nFaces
    iPatch = pC2f(1,ifl)
    ifg    = pC2f(2,ifl)

! DEBUG
!    WRITE(0,*) '@@@200',pRegion%iRegionGlobal,'Comp xsect with face:',iPatch,ifg
!        IF ( iPatch == 0 ) THEN
!          WRITE(0,*) '@@@210',pRegion%iRegionGlobal,pGrid%fn(1:3,ifg)
!        ELSE
!          pPatch => pRegion%patches(iPatch)
!          WRITE(0,*) '@@@211',pPatch%fn(1:3,ifg)         
!        END IF ! iPatch
! END DEBUG

! ==============================================================================
!   Interior face
! ==============================================================================

    IF ( iPatch == 0 ) THEN 
    
! ------------------------------------------------------------------------------
!     Triangular face
! ------------------------------------------------------------------------------    
    
! TEMPORARY
!      IF ( pGrid%f2v(4,ifg) == VERT_NONE ) THEN 
! END TEMPORARY
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

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

        xCofg = pGrid%fc(XCOORD,ifg)
        yCofg = pGrid%fc(YCOORD,ifg)
        zCofg = pGrid%fc(ZCOORD,ifg) 

        numer(ifl) = (xCofg-xLoc)*fnx + (yCofg-yLoc)*fny + (zCofg-zLoc)*fnz
        denom(ifl) = ex*fnx + ey*fny + ez*fnz 

! ------------------------------------------------------------------------------
!     Quadrilateral face
! ------------------------------------------------------------------------------    
    
! TEMPORARY
!      ELSE 
!        CALL RFLU_BLIN_ComputeXSectLine(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
!                                        iPatch,ifg,nXSect,xsd)        
!
!! DEBUG
!        WRITE(0,*) '@@@300',pRegion%iRegionGlobal,'xsect dist:',nXSect,xsd
!! END DEBUG
!                              
!        IF ( nXSect > 0 ) THEN 
!          numer(ifl) = xsd(1)
!          denom(ifl) = 1.0_RFREAL
!        ELSE 
!          numer(ifl) = HUGE(1.0_RFREAL)
!          denom(ifl) = 1.0_RFREAL 
!        END IF ! nXSect  
!      END IF ! pGrid%f2v     
! END TEMPORARY

! ==============================================================================
!   Boundary face
! ==============================================================================

    ELSE IF ( iPatch > 0 ) THEN 
      pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!     Triangular face
! ------------------------------------------------------------------------------    
    
!      IF ( pPatch%bf2v(4,ifg) == VERT_NONE ) THEN 
        fnx = pPatch%fn(XCOORD,ifg)
        fny = pPatch%fn(YCOORD,ifg)
        fnz = pPatch%fn(ZCOORD,ifg)          

        xCofg = pPatch%fc(XCOORD,ifg)
        yCofg = pPatch%fc(YCOORD,ifg)
        zCofg = pPatch%fc(ZCOORD,ifg)
      
        numer(ifl) = (xCofg-xLoc)*fnx + (yCofg-yLoc)*fny + (zCofg-zLoc)*fnz
        denom(ifl) = ex*fnx + ey*fny + ez*fnz      

! ------------------------------------------------------------------------------
!     Quadrilateral face
! ------------------------------------------------------------------------------    
    
! TEMPORARY
!      ELSE 
!        CALL RFLU_BLIN_ComputeXSectLine(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
!                                        iPatch,ifg,nXSect,xsd)        
!                              
!        IF ( nXSect > 0 ) THEN 
!          numer(ifl) = xsd(1)
!          denom(ifl) = 1.0_RFREAL
!        ELSE 
!          numer(ifl) = HUGE(1.0_RFREAL)
!          denom(ifl) = 1.0_RFREAL 
!        END IF ! nXSect      
!      END IF ! pPatch%bf2v
! END TEMPORARY
    ELSE ! Defensive programming
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iPatch
  END DO ! ifl

! ******************************************************************************
! Set output
! ******************************************************************************

  errorFlag = ERR_NONE

  distOut1 = HUGE(1.0_RFREAL)
  distOut2 = HUGE(1.0_RFREAL)

  iflOut1 = CRAZY_VALUE_INT
  iflOut2 = CRAZY_VALUE_INT

  faceLoop: DO ifl = 1,nFaces
    IF ( numer(ifl) > toler ) THEN ! Inside cell icl
      IF ( denom(ifl) > 0.0 ) THEN
        distOut = MAX(numer(ifl),0.0_RFREAL)/denom(ifl)
        
        IF ( distOut < distOut1 ) THEN 
          distOut2 = distOut1        
          iflOut2  = iflOut1
          
          distOut1 = distOut
          iflOut1  = ifl
        ELSE IF ( distOut < distOut2 ) THEN 
          distOut2 = distOut
          iflOut2  = ifl        
        END IF ! dist
      END IF ! denom 
    ELSE ! Outside cell icl
      errorFlag = 1
      distOut = 0.0_RFREAL
      iflOut  = ifl
      
      EXIT faceLoop
    END IF ! numer
  END DO faceLoop

  IF ( errorFlag == ERR_NONE ) THEN
    distOut = distOut1 ! Temporary 
    iflOut  = iflOut1     
  END IF ! errorFlag

! DEBUG
!  WRITE(0,*) '@@@400',pRegion%iRegionGlobal,iflOut,distOut
! END DEBUG
                              
  iPatchOut = pC2f(1,iflOut)
  ifgOut    = pC2f(2,iflOut)      

  xLoc = xLoc + distOut*ex
  yLoc = yLoc + distOut*ey
  zLoc = zLoc + distOut*ez    

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeLineCellXSectSafe






! ******************************************************************************
!
! Purpose: Determine whether point is inside bounding box.
!
! Description: 
!
! Input:
!   global              Pointer to global data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   xMin,yMin,zMin      Minimum x-, y-, and z-coordinates
!   xMax,yMax,zMax      Maximum x-, y-, and z-coordinates
!
! Output: 
!   RFLU_TestInBoundingBox = .TRUE.     if point inside bounding box
!   RFLU_TestInBoundingBox = .FALSE.    if point outside bounding box
!
! Notes: None.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax, & 
                                     yMin,yMax,zMin,zMax)
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: xLoc,xMax,xMin,yLoc,yMax,yMin,zLoc,zMax,zMin
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TestInBoundBox',&
  'RFLU_ModGeometryTools.F90')

! ******************************************************************************
! Determine whether point is in bounding box
! ******************************************************************************

  RFLU_TestInBoundBox = .TRUE. 

  IF ( (xLoc < xMin) .OR. (xLoc > xMax) ) THEN 
    RFLU_TestInBoundBox = .FALSE. 
    CALL DeregisterFunction(global)    
    RETURN 
  END IF ! xLoc

  IF ( (yLoc < yMin) .OR. (yLoc > yMax) ) THEN 
    RFLU_TestInBoundBox = .FALSE. 
    CALL DeregisterFunction(global)      
    RETURN 
  END IF ! yLoc

  IF ( (zLoc < zMin) .OR. (zLoc > zMax) ) THEN 
    RFLU_TestInBoundBox = .FALSE. 
    CALL DeregisterFunction(global)      
    RETURN 
  END IF ! zLoc

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_TestInBoundBox





! ******************************************************************************
!
! Purpose: Determine whether vector is aligned with Cartesian coordinate 
!   direction.
!
! Description: 
!
! Input:
!   global              Pointer to global data
!   dr			Vector
!   dir			Cartesian direction index
!
! Output: 
!   RFLU_TestVectorCartAxisAligned = .TRUE.     if dr aligned with Cartesian 
!                                               coordinate direction
!   RFLU_TestVectorCartAxisAligned = .FALSE.    if dr not aligned with Cartesian 
!                                      		coordinate direction
!
! Notes: None.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_TestVectorCartAxisAligned(global,dr,dir)
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: dir
  REAL(RFREAL), INTENT(INOUT) :: dr(XCOORD:ZCOORD)
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: drSum,drTol

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_TestVectorCartAxisAligned',&
  'RFLU_ModGeometryTools.F90')

! ******************************************************************************
! Set tolerance
! ******************************************************************************

  drTol = 1.0E-12_RFREAL

! ******************************************************************************
! Test whether vector aligned 
! ******************************************************************************

  RFLU_TestVectorCartAxisAligned = .FALSE. 

  dr(dir) = 0.0_RFREAL
        
  drSum = ABS(dr(XCOORD)) + ABS(dr(YCOORD)) + ABS(dr(ZCOORD))

  IF ( drSum <= drTol ) THEN 
    RFLU_TestVectorCartAxisAligned = .TRUE.
  END IF ! drSum

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_TestVectorCartAxisAligned





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModGeometryTools


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGeometryTools.F90,v $
! Revision 1.5  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 18:06:42  haselbac
! Added function to test for alignment ofvector with Cartesian direction
!
! Revision 1.2  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.1  2005/12/24 21:17:50  haselbac
! Initial revision
!
! ******************************************************************************
  










