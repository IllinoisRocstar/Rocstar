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
! Purpose: Suite of routines for particle tracking on Eulerian grid.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_ModFindCells.F90,v 1.16 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_RFLU_ModFindCells
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModPartLag,    ONLY: t_plag, t_plag_input
  USE ModBndPatch,   ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModBorder,     ONLY: t_border
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters
    
  USE ModTools, ONLY: FloatEqual
  
  USE RFLU_ModCellFaceEdgeInfo, ONLY: RFLU_GetGlobalCellKind, & 
                                      RFLU_GetFaceKind  
  USE RFLU_ModGeometryTools
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell, & 
                                RFLU_ICT_TestInCellFancy, & 
                                RFLU_ICT_TestInCellLohner
  USE RFLU_ModOctree
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformVector

  USE PLAG_ModInterfaces, ONLY: PLAG_ReflectParticleData 
  USE PLAG_ModSurfStats, ONLY: PLAG_GatherSurfStats

  USE RFLU_ModMPI, ONLY: RFLU_MPI_RecreateBufferIPclSend  

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_RFLU_ComputeDistTot,    &
            PLAG_RFLU_FindCellsBrute,    &
            PLAG_RFLU_FindCellsBruteMod, &
            PLAG_RFLU_FindCellsLohner,   & 
            PLAG_RFLU_FindCellsOct,      &
            PLAG_RFLU_FindCellsOctMod,   &
            PLAG_RFLU_FindCellsTrajFast, &
            PLAG_RFLU_FindCellsTrajSafe
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_ModFindCells.F90,v $ $Revision: 1.16 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS











! ******************************************************************************
!
! Purpose: Compute total distance travelled.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeDistTot(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: errorString
  INTEGER ::iPcl
  REAL(RFREAL) :: xLocNew,xLocOld,xTraj,yLocNew,yLocOld,yTraj,zLocNew, &
                  zLocOld,zTraj
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ComputeDistTot',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls
 
! ==============================================================================  
!   Compute distance travelled and trajectory
! ==============================================================================  
  
    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
        
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld      
    zTraj = zLocNew - zLocOld

    pPlag%arv(ARV_PLAG_DISTOT,iPcl)  = SQRT( xTraj*xTraj + yTraj*yTraj &
                                           + zTraj*zTraj )
          
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_ComputeDistTot









! ******************************************************************************
!
! Purpose: Determine cells which contain particles using brute force approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBrute(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: errorFlag,icg,iPcl
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBrute',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls  

! ==============================================================================  
!   Get particle position
! ==============================================================================  

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Loop over actual cells
! ==============================================================================  
  
    foundFlag = .FALSE.
  
    cellLoop: DO icg = 1,pGrid%nCells         
      IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
        pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        foundFlag = .TRUE.

        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell        
    END DO cellLoop
      
! ==============================================================================  
!   Check whether particle was found
! ==============================================================================  
            
    IF ( foundFlag .EQV. .FALSE. ) THEN
      WRITE(errorString,'(I6)') iPcl 
      CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))
    END IF ! foundFlag
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBrute







! ******************************************************************************
!
! Purpose: Kernel for brute force search.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   xLoc,yLoc,zLoc      Location
!
! Output:
!   icgOut              Cell containing location
!
! Notes: 
!   1. icgOut is set to CRAZY_VALUE_INT if no cell contained location.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBruteKernel(pRegion,xLoc,yLoc,zLoc,icgOut)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(OUT) :: icgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag,icg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBruteKernel',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Loop over cells
! ******************************************************************************

  icgOut = CRAZY_VALUE_INT

  cellLoop: DO icg = 1,pGrid%nCells         
    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
      icgOut = icg

      EXIT cellLoop
    END IF ! RFLU_ICT_TestInCell        
  END DO cellLoop
      
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBruteKernel







! ******************************************************************************
!
! Purpose: Determine cells which contain particles using modified brute 
!  force approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   2. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBruteMod(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: errorString
  INTEGER :: errorFlag,icg,iPcl
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBruteMod',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls  

! ==============================================================================  
!   Get particle cell and position
! ==============================================================================  

    icg  = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl) ! NOTE OLD cell index

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl) ! NOTE already NEW particle position
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Check if particle still in same cell 
! ==============================================================================  

    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .FALSE. ) THEN    

! ------------------------------------------------------------------------------  
!     Find cell containing particle
! ------------------------------------------------------------------------------  
  
      CALL PLAG_RFLU_FindCellsBruteKernel(pRegion,xLoc,yLoc,zLoc,icg)
  
      IF ( icg /= CRAZY_VALUE_INT ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
      ELSE 

! TEMPORARY
        WRITE(*,*) 'timeCurrent,iPcl,xLocOld,yLocOld,zLocOld,xLoc,yLoc,zLoc = ',&
        global%currentTime,iPcl,pPlag%cvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcl),&
        xLoc,yLoc,zLoc 
! END TEMPORARY

        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))     
      END IF ! icg
            
! ==============================================================================  
!   Particle still in same cell, copy cell index
! ==============================================================================        
      
    ELSE   
      pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    END IF ! RFLU_ICT_TestInCell
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBruteMod







! ******************************************************************************
!
! Purpose: Determine cells which contain particles using Lohners approach.
!
! Description: Follow particle path by passing particle to face which failed
!   in-cell test by largest amount.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. See R. Lohner and J. Ambrosiano, A Vectorized Particle Tracer for 
!      Unstructured Grids, J. Comp. Phys., Vol. 91, pp. 22-31, 1990.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsLohner(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: testInCell
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,icg,icgOut,ifg,iloc,iPcl,loopCounter
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsLohner',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls
    loopCounter = 0 ! Reset loop counter 

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLoc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)  
                   
    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    infLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
                    
      CALL RFLU_ICT_TestInCellLohner(pRegion,xLoc,yLoc,zLoc,icg,testInCell, &
                                     iloc,ifg)                                                        
                                   
! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- In-cell test successful --------------------------------------------------

      IF ( testInCell .EQV. .TRUE. ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        EXIT infLoop
        
! --- In-cell test failed ------------------------------------------------------
        
      ELSE 

! ----- Intersect interior face

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1 
! TO DO
!            CASE ( FACE_KIND_AV ) ! Actual-virtual face
! Count particle intersections for each buffer           
! END TO DO
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Intersect boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)

          SELECT CASE ( pPatch%bcType )             
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType
        END IF ! iLoc 
      END IF ! testInCell
      
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
                    
    END DO infLoop     
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsLohner







! ******************************************************************************
!
! Purpose: Determine cells which contain particles using Octree approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cells which contain an 
!      initial distribution of particles, but should not be used to track 
!      moving particles, because it will be inefficient on account of it
!      not using knowledge of past position. It will also overwrite entry 
!      containing initial region in particle data structure. 
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOct(pRegion)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: ccSize,errorFlag,icg,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax,yMin,zDel,zLoc, & 
                  zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOct',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  delFrac = 0.01_RFREAL
  
  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************  
! Build bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel 
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel 
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel  

! ******************************************************************************  
! Build Octree with cell centroids: NOTE only for actual cells
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)

  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls

! ==============================================================================  
!   Test particle location against bounding box of partition
! ==============================================================================  
  
    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)
  
    IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                             zMin,zMax) .EQV. .TRUE. ) THEN  
                               
! ------------------------------------------------------------------------------  
!     Query octree to get closest cells  
! ------------------------------------------------------------------------------ 
    
      CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)
        
! ------------------------------------------------------------------------------
!     Test cells obtained from Octree if they contain specified location
! ------------------------------------------------------------------------------        
           
      foundFlag = .FALSE.     
           
      cellLoop: DO icg = 1,ccSize           
        IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icg)) .EQV. .TRUE. ) THEN
          pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
          pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = cc(icg)

          foundFlag = .TRUE.

          EXIT cellLoop
        END IF ! RFLU_ICT_TestInCell        
      END DO cellLoop
      
! ------------------------------------------------------------------------------
!     Check whether particle was found
! ------------------------------------------------------------------------------        
            
      IF ( foundFlag .EQV. .FALSE. ) THEN
        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, &
                       TRIM(errorString))
      END IF ! foundFlag

! ==============================================================================  
!   Particle located outside bounding box of partition
! ==============================================================================  

    ELSE   
      WRITE(errorString,'(I6)') iPcl 
      CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))      
    END IF ! RFLU_TestInBoundNox 
  END DO ! iPcl

! ******************************************************************************
! Destroy Octree and deallocate temporary memory
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOct








! ******************************************************************************
!
! Purpose: Kernel for Octree search
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   xLoc,yLoc,zLoc      Location
!
! Output:
!   icgOut              Cell containing location
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   3. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOctKernel(pRegion,xLoc,yLoc,zLoc,xMin,xMax,&
                                        yMin,yMax,zMin,zMax,icgOut)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(INOUT)   :: icgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,xMax,xMin,yLoc,yMax,yMin,zLoc,zMax,zMin
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  INTEGER :: ccSize,errorFlag,icgLoop,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOctKernel',&
  'PLAG_RFLU_ModFindCells.F90')
  
! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid

!  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

! TEMPORARY
  ccSize = MIN(50,pGrid%nCells) ! Must be larger than unity 
! END TEMPORARY

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! Test particle location against bounding box of partition
! ******************************************************************************

  IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                             zMin,zMax) .EQV. .TRUE. ) THEN  

! ============================================================================== 
!  Query octree to get closest cells  
! ============================================================================== 
    
   CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)

! ============================================================================== 
!  Test cells obtained from Octree if they contain specified location
! ============================================================================== 
! ------------------------------------------------------------------------------ 
           
   foundFlag = .FALSE.     

   cellLoop: DO icgLoop = 1,ccSize           
      IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icgLoop)) .EQV. .TRUE. ) THEN
        icgOut = cc(icgLoop)

        foundFlag = .TRUE.

        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell        
    END DO cellLoop

! ============================================================================== 
!   Check whether particle was found 
! ============================================================================== 
            
    IF ( foundFlag .EQV. .FALSE. ) THEN
      icgOut = CRAZY_VALUE_INT
    END IF ! foundFlag

! ******************************************************************************
! Particle located outside bounding box of partition
! ******************************************************************************

  ELSE   
    icgOut = CRAZY_VALUE_INT    
  END IF ! RFLU_TestInBoundNox

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOctKernel








! ******************************************************************************
!
! Purpose: Determine cells which contain particles using modified Octree 
!   approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   3. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOctMod(pRegion)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: ccSize,errorFlag,icg,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax,yMin,zDel,zLoc, & 
                  zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOctMod',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  delFrac = 0.01_RFREAL
  
!  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

! TEMPORARY
  ccSize = MIN(50,pGrid%nCells) ! Must be larger than unity 
! END TEMPORARY

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************  
! Build bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel 
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel 
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel  

! ******************************************************************************  
! Build Octree with cell centroids: NOTE only for actual cells
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)

  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls

! ==============================================================================  
!   Get particle cell and position
! ==============================================================================  

    icg  = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl) ! NOTE OLD cell index
  
    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl) ! NOTE already NEW particle position
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Check if particle still in same cell 
! ==============================================================================  
    
    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .FALSE. ) THEN    
 
! ------------------------------------------------------------------------------  
!     Test particle location against bounding box of partition
! ------------------------------------------------------------------------------  
  
      IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                               zMin,zMax) .EQV. .TRUE. ) THEN  
                               
! ----- Query octree to get closest cells -------------------------------------- 
    
        CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)
        
! ----- Test cells obtained from Octree if they contain specified location -----
           
        foundFlag = .FALSE.     

        cellLoop: DO icg = 1,ccSize           
          IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icg)) .EQV. .TRUE. ) THEN
            pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = cc(icg)

            foundFlag = .TRUE.

            EXIT cellLoop
          END IF ! RFLU_ICT_TestInCell        
        END DO cellLoop
      
! ----- Check whether particle was found ---------------------------------------
            
        IF ( foundFlag .EQV. .FALSE. ) THEN
          WRITE(errorString,'(I6)') iPcl 
          CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, &
                         TRIM(errorString))
        END IF ! foundFlag

! ------------------------------------------------------------------------------  
!     Particle located outside bounding box of partition
! ------------------------------------------------------------------------------  

      ELSE   

! TEMPORARY
        WRITE(*,*) 'timeCurrent,iPcl,xLocOld,yLocOld,zLocOld,xLoc,yLoc,zLoc = ',&
        global%currentTime,iPcl,pPlag%cvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcl),&
        xLoc,yLoc,zLoc
! END TEMPORARY

        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))      
      END IF ! RFLU_TestInBoundNox

! ==============================================================================  
!   Particle still in same cell, copy cell index 
! ==============================================================================  
    
    ELSE   
      pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    END IF ! RFLU_ICT_TestInCell 
  END DO ! iPcl

! ******************************************************************************
! Destroy Octree and deallocate temporary memory
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOctMod







! ******************************************************************************
!
! Purpose: Determine cells which contain particles by following trajectory 
!   using fast algorithm.
!
! Description: Follow particle path by intersecting particle trajectory with
!   faces and taking appropriate action determined by type of intersected face.
!
! Input: 
!  pRegion      Pointer to region
!  iPclBeg      Beginning particle index
!  iPclEnd      Ending particle index
!
! Output: None.
!
! Notes:
!   1. Need to properly handle particle deletion from outflow or motion
!      between domain during an intermediate RK-stage as its final 
!      position at (n+1) step might not be outside the region.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsTrajFast(pRegion,iPclBeg,iPclEnd)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER :: iPclBeg,iPclEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: inCellCheckFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,iBorder,icg,ifg,ifgOut,ifl,iloc,ilocOut, & 
             iPatch,iPcl,loopCounter
  REAL(RFREAL) :: dist,distTot,distTotCutoff,eps,fnx,fny,fnz,iMagTraj,theta, &
                  xLoc,xLocNew,xLocOld,xTraj,yLoc,yLocNew,yLocOld,yTraj,zLoc, &
                  zLocNew,zLocOld,zTraj
  TYPE(t_border), POINTER :: pBorder
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsTrajFast',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  eps = EPSILON(1.0_RFREAL)
  distTotCutoff = 10*EPSILON(1.0_RFREAL) ! NOTE must be less than tolerICT

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = iPclBeg,iPclEnd
    loopCounter = 0 ! Reset loop counter 

! ==============================================================================  
!   Set distance travelled and compute trajectory. NOTE that cannot compute 
!   trajectory (unit) vector from distTot because distTot diminishes as the 
!   particle travels along the trajectory but xLocNew and xLocOld (as well as
!   y and z components) do not change.
! ==============================================================================  

    distTot = pPlag%arv(ARV_PLAG_DISTOT,iPcl)

    IF ( distTot < distTotCutOff ) THEN  
      CYCLE 
    END IF ! distTot

    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
        
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld     
    zTraj = zLocNew - zLocOld
             
    iMagTraj = 1.0_RFREAL/(SQRT(xTraj*xTraj + yTraj*yTraj + zTraj*zTraj) + eps)
    
    xTraj = iMagTraj*xTraj
    yTraj = iMagTraj*yTraj
    zTraj = iMagTraj*zTraj    
    
! ------------------------------------------------------------------------------
!   Set variables for trajectory loop
! ------------------------------------------------------------------------------    
    
    xLoc = xLocNew - distTot*xTraj
    yLoc = yLocNew - distTot*yTraj
    zLoc = zLocNew - distTot*zTraj

    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    trajLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
                    
      CALL RFLU_ComputeLineCellXSectFast(pRegion,xLoc,yLoc,zLoc,xTraj,yTraj, &
                                         zTraj,icg,dist,iloc,ifg)

! ------------------------------------------------------------------------------
!     Update total distance travelled                             
! ------------------------------------------------------------------------------
                                   
      distTot = distTot - dist

! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- No distance remaining: Found cell containing new location ----------------

      IF ( distTot <= 0.0_RFREAL ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        EXIT trajLoop
        
! --- Distance remaining: Keep searching ---------------------------------------        
        
      ELSE 

! ----- Trajectory intersects interior face

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1                                    
            CASE ( FACE_KIND_AV ) ! Actual-virtual face
              ifl = ifg - pGrid%nFaces + pGrid%nFacesAV

              iBorder = pGrid%avf2b(1,ifl)

              pBorder => pGrid%borders(iBorder)
              
              pBorder%nPclsSend = pBorder%nPclsSend + 1

              IF ( pBorder%nPclsSend > pBorder%nPclsSendMax ) &
                CALL RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)
          
              pBorder%iPclSend(1,pBorder%nPclsSend) = iPcl
              pBorder%iPclSend(2,pBorder%nPclsSend) = pGrid%avf2b(2,ifl)            

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_COMM
              pPlag%arv(ARV_PLAG_DISTOT,iPcl) = distTot

              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1
              
              pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
              
	      iPatch = pGrid%avf2p(ifl)
	      
	      IF ( iPatch /= CRAZY_VALUE_INT ) THEN ! Transform coordinates
	        pPatch => pRegion%patches(iPatch)
		
		IF ( pPatch%bcType == BC_PERIODIC ) THEN 
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocOld, &
		                                 yLocOld,zLocOld)
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocNew, &
		                                 yLocNew,zLocNew)						 
		ELSE
		  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
		END IF ! pPatch%bcType
	      END IF ! iPatch	      
	      
              EXIT trajLoop                  
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Trajectory intersects boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)
          
          fnx = pPatch%fn(XCOORD,ifg)
          fny = pPatch%fn(YCOORD,ifg)
          fnz = pPatch%fn(ZCOORD,ifg) 
                      
          theta = xTraj*fnx + yTraj*fny + zTraj*fnz 
          
          SELECT CASE ( pPatch%bcType )             
            CASE ( BC_SLIPWALL:BC_SLIPWALL+BC_RANGE )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_NOSLIPWALL:BC_NOSLIPWALL+BC_RANGE ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_INJECTION:BC_INJECTION+BC_RANGE ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_OUTFLOW:BC_OUTFLOW+BC_RANGE )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop            
            CASE ( BC_FARFIELD:BC_FARFIELD+BC_RANGE )
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop 
            CASE ( BC_SYMMETRY:BC_SYMMETRY+BC_RANGE )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_VIRTUAL:BC_VIRTUAL+BC_RANGE )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType

        END IF ! iLoc 
      END IF ! distTot
    
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
             'Infinite loop encountered in particle cell search algorithm'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).'

        WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime

        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                         '    iPcl     ', &
                                         '    idIni    ', &
                                         '    RegIni   ', &
                                         '    icg      ', &
                                         '  x-location ', &
                                         '  y-location ', &
                                         '  z-Location ', &
                                         '   Energy    ', &
                                         '   Diameter  '

        WRITE(STDOUT,'(A,4X,4(1X,I8),6(1X,E13.6))') SOLVER_NAME,iPcl, &
                        pPlag%aivOld(AIV_PLAG_PIDINI,iPcl),           &
                        pPlag%aivOld(AIV_PLAG_REGINI,iPcl),           &
                        icg,xLoc,yLoc,zLoc,                           &
                        pPlag%cv(CV_PLAG_ENER,iPcl),                  &
                        pPlag%dv(DV_PLAG_DIAM,iPcl)

        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter        
    END DO trajLoop 

! ==============================================================================
!   Check that kept particles are indeed located in new cells
! ==============================================================================        

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
	 (pPlag%aiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_KEEP) ) THEN
      CALL RFLU_ICT_TestInCellFancy(pRegion,xLocNew,yLocNew,zLocNew,icg, &
				    inCellCheckFlag,ilocOut,ifgOut)

      IF ( inCellCheckFlag .EQV. .FALSE. ) THEN 
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Particle index:',iPcl
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Cell which failed test:', & 
						   icg
	WRITE(STDERR,'(A,1X,A,2(1X,I6))') SOLVER_NAME, &
					  'Face which failed test:', & 
					  ilocOut,ifgOut								 
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle old location:', &
					      xLocOld,yLocOld,zLocOld
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle new location:', &
					      xLocNew,yLocNew,zLocNew

	WRITE(errorString,'(I6)') iPcl
	CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, & 
		       TRIM(errorString))	    
      END IF ! inCellCheckFlag
    END IF ! global%checkLevel    
    
! ==============================================================================
!   Store new position. IMPORTANT: Need to store new position because it may 
!   have been reflected and old position because it may have been transformed.
! ==============================================================================

    pPlag%cv(CV_PLAG_XPOS,iPcl) = xLocNew
    pPlag%cv(CV_PLAG_YPOS,iPcl) = yLocNew
    pPlag%cv(CV_PLAG_ZPOS,iPcl) = zLocNew
     
    pPlag%cvOld(CV_PLAG_XPOS,iPcl) = xLocOld
    pPlag%cvOld(CV_PLAG_YPOS,iPcl) = yLocOld
    pPlag%cvOld(CV_PLAG_ZPOS,iPcl) = zLocOld       
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsTrajFast







! ******************************************************************************
!
! Purpose: Determine cells which contain particles by following trajectory 
!   using safe algorithm.
!
! Description: Follow particle path by intersecting particle trajectory with
!   faces and taking appropriate action determined by type of intersected face.
!
! Input: 
!  pRegion      Pointer to region
!  iPclBeg      Beginning particle index
!  iPclEnd      Ending particle index
!
! Output: None.
!
! Notes:
!   1. Need to properly handle particle deletion from outflow or motion
!      between domain during an intermediate RK-stage as its final 
!      position at (n+1) step might not be outside the region.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsTrajSafe(pRegion,iPclBeg,iPclEnd)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER :: iPclBeg,iPclEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: inCellCheckFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,iBorder,icg,ifg,ifgOut,ifl,iloc, &
             ilocOut,iPatch,iPcl,loopCounter
  REAL(RFREAL) :: dist,distTot,distTotCutoff,eps,fnx,fny,fnz,iMagTraj,theta, &
                  xLoc,xLocNew,xLocOld,xTraj,yLoc,yLocNew,yLocOld,yTraj,zLoc, &
                  zLocNew,zLocOld,zTraj
  TYPE(t_border), POINTER :: pBorder
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsTrajSafe',&
  'PLAG_RFLU_ModFindCells.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  eps = EPSILON(1.0_RFREAL)
  distTotCutoff = 10*eps ! NOTE must be less than tolerICT

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = iPclBeg,iPclEnd
    loopCounter = 0 ! Reset loop counter 

! ==============================================================================  
!   Set distance travelled and compute trajectory. NOTE that cannot compute 
!   trajectory (unit) vector from distTot because distTot diminishes as the 
!   particle travels along the trajectory but xLocNew and xLocOld (as well as
!   y and z components) do not change.
! ==============================================================================  

    distTot = pPlag%arv(ARV_PLAG_DISTOT,iPcl)
    
    IF ( distTot < distTotCutoff ) THEN  
      CYCLE 
    END IF ! distTot

    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
        
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld      
    zTraj = zLocNew - zLocOld

    iMagTraj = 1.0_RFREAL/(SQRT(xTraj*xTraj + yTraj*yTraj + zTraj*zTraj) + eps)
    
    xTraj = iMagTraj*xTraj
    yTraj = iMagTraj*yTraj
    zTraj = iMagTraj*zTraj    
    
! ------------------------------------------------------------------------------
!   Set variables for trajectory loop
! ------------------------------------------------------------------------------    
    
    xLoc = xLocNew - distTot*xTraj
    yLoc = yLocNew - distTot*yTraj
    zLoc = zLocNew - distTot*zTraj
    
    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    trajLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
                    
      CALL RFLU_ComputeLineCellXSectSafe(pRegion,xLoc,yLoc,zLoc,xTraj,yTraj, &
                                         zTraj,icg,dist,iloc,ifg)  
                                                                                              
! ------------------------------------------------------------------------------
!     Update total distance travelled                             
! ------------------------------------------------------------------------------
                                   
      distTot = distTot - dist
      
! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- No distance remaining: Found cell containing new location ----------------

      IF ( distTot <= 0.0_RFREAL ) THEN       
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
	
        EXIT trajLoop
        
! --- Distance remaining: Keep searching ---------------------------------------        
        
      ELSE 

! ----- Trajectory intersects interior face

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1         
            CASE ( FACE_KIND_AV ) ! Actual-virtual face
              ifl = ifg - pGrid%nFaces + pGrid%nFacesAV

              iBorder = pGrid%avf2b(1,ifl)

              pBorder => pGrid%borders(iBorder)
              
              pBorder%nPclsSend = pBorder%nPclsSend + 1

              IF ( pBorder%nPclsSend > pBorder%nPclsSendMax ) &
                CALL RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)
              
              pBorder%iPclSend(1,pBorder%nPclsSend) = iPcl
              pBorder%iPclSend(2,pBorder%nPclsSend) = pGrid%avf2b(2,ifl)            

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_COMM
              pPlag%arv(ARV_PLAG_DISTOT,iPcl) = distTot
              
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1
              
              pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg  
              
	      iPatch = pGrid%avf2p(ifl)
	      
	      IF ( iPatch /= CRAZY_VALUE_INT ) THEN ! Transform coordinates
	        pPatch => pRegion%patches(iPatch)
		
		IF ( pPatch%bcType == BC_PERIODIC ) THEN 
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocOld, &
		                                 yLocOld,zLocOld)
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocNew, &
		                                 yLocNew,zLocNew)						 
		ELSE
		  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
		END IF ! pPatch%bcType
	      END IF ! iPatch

              EXIT trajLoop                             
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Trajectory intersects boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)
          
          fnx = pPatch%fn(XCOORD,ifg)
          fny = pPatch%fn(YCOORD,ifg)
          fnz = pPatch%fn(ZCOORD,ifg) 
                      
          theta = xTraj*fnx + yTraj*fny + zTraj*fnz 
          
          SELECT CASE ( pPatch%bcType )             
            CASE ( BC_SLIPWALL:BC_SLIPWALL+BC_RANGE )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_NOSLIPWALL:BC_NOSLIPWALL+BC_RANGE ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_INJECTION:BC_INJECTION+BC_RANGE ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_OUTFLOW:BC_OUTFLOW+BC_RANGE )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop            
            CASE ( BC_FARFIELD:BC_FARFIELD+BC_RANGE )
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop 
            CASE ( BC_SYMMETRY:BC_SYMMETRY+BC_RANGE )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_VIRTUAL:BC_VIRTUAL+BC_RANGE )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType

        END IF ! iLoc 
      END IF ! distTot
    
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
             'Infinite loop encountered in particle cell search algorithm'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).'        
          
        WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime              
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
        WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                         '    iPcl     ', &
                                         '    idIni    ', &
                                         '    RegIni   ', &
                                         '    icg      ', &
                                         '  x-location ', &
                                         '  y-location ', &
                                         '  z-Location ', &
                                         '   Energy    ', &
                                         '   Diameter  '       

        WRITE(STDOUT,'(A,4X,4(1X,I8),6(1X,E13.6))') SOLVER_NAME,iPcl, & 
                        pPlag%aivOld(AIV_PLAG_PIDINI,iPcl),           &
                        pPlag%aivOld(AIV_PLAG_REGINI,iPcl),           &
                        icg,xLoc,yLoc,zLoc,                           &
                        pPlag%cv(CV_PLAG_ENER,iPcl),                  &
                        pPlag%dv(DV_PLAG_DIAM,iPcl)

        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter        
    END DO trajLoop 

! ==============================================================================
!   Check that kept particles are indeed located in new cells
! ==============================================================================        

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
	 (pPlag%aiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_KEEP) ) THEN
      CALL RFLU_ICT_TestInCellFancy(pRegion,xLocNew,yLocNew,zLocNew,icg, &
				    inCellCheckFlag,ilocOut,ifgOut)

      IF ( inCellCheckFlag .EQV. .FALSE. ) THEN 
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Particle index:',iPcl
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Cell which failed test:', & 
						   icg
	WRITE(STDERR,'(A,1X,A,2(1X,I6))') SOLVER_NAME, &
					  'Face which failed test:', & 
					  ilocOut,ifgOut								 
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle old location:', &
					      xLocOld,yLocOld,zLocOld
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle new location:', &
					      xLocNew,yLocNew,zLocNew

	WRITE(errorString,'(I6)') iPcl
	CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, & 
		       TRIM(errorString))	    
      END IF ! inCellCheckFlag
    END IF ! global%checkLevel    
    
! ==============================================================================
!   Store new position. IMPORTANT: Need to store new position because it may 
!   have been reflected and old position because it may have been transformed.
! ==============================================================================

    pPlag%cv(CV_PLAG_XPOS,iPcl) = xLocNew
    pPlag%cv(CV_PLAG_YPOS,iPcl) = yLocNew
    pPlag%cv(CV_PLAG_ZPOS,iPcl) = zLocNew
     
    pPlag%cvOld(CV_PLAG_XPOS,iPcl) = xLocOld
    pPlag%cvOld(CV_PLAG_YPOS,iPcl) = yLocOld
    pPlag%cvOld(CV_PLAG_ZPOS,iPcl) = zLocOld                  	      
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsTrajSafe





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_RFLU_ModFindCells

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ModFindCells.F90,v $
! Revision 1.16  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2006/08/18 21:11:35  fnajjar
! Enabled periodicity, cosmetics
!
! Revision 1.13  2006/05/05 02:15:49  haselbac
! Bug fixes: Incorrect computation of traj, missing EXIT, wrong IF
!
! Revision 1.12  2006/05/02 17:46:23  fnajjar
! Allowed surface statistics to be gathered on active patches
!
! Revision 1.11  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.10  2005/12/24 21:39:47  haselbac
! Adapted to changes in ICT
!
! Revision 1.9  2005/12/19 16:47:45  fnajjar
! Added verbosity around error trap for infinite loop counter
!
! Revision 1.8  2005/12/14 21:21:21  fnajjar
! Added call for dynamic allocation of iPclSend
!
! Revision 1.7  2005/12/02 20:07:50  fnajjar
! Added particle reflection for BC_VIRTUAL
!
! Revision 1.6  2005/11/12 00:34:08  fnajjar
! Bug fix for iMagTraj in Fast and Safe
!
! Revision 1.5  2005/09/20 15:46:35  fnajjar
! Fixed bug in TrajFast and added error trap for iPclSend overflow
!
! Revision 1.4  2005/07/18 20:47:41  fnajjar
! Aligned fast routine with safe adding proper construct for trajectory
!
! Revision 1.3  2005/05/18 22:17:50  fnajjar
! Added PLAG_RFLU_ComputeDistTot, changed comp of distTot and xyzLoc, 
! infrast for comm
!
! Revision 1.2  2005/05/02 21:59:38  haselbac
! Preparation for parallelization of FindCells routines, cosmetics
!
! Revision 1.1  2005/04/27 14:55:55  fnajjar
! Initial import
!
! ******************************************************************************
















