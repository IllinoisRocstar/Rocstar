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
! Purpose: Set options for moving grid computations.
!
! Description: The options for imposing boundary conditions on the grid motion
!   are determined by whether the patches are flat and whether they are aligned
!   with a coordinate direction. If not flat, no boundary condition is (can be)
!   applied without a geometrical description of the surface. If the patch is
!   flat, Dirichlet conditions are set if the patch normal is perfectly aligned
!   with a coordinate direction, otherwise von Neumann conditions are applied.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: 
!   1. At present, the only option which is determined in this routine is 
!      how the boundary conditions are set for the grid motion, namely either
!      von Neumann or Dirichlet conditions. The von Neumann condition is 
!      set for those patches the normal to which is not aligned with a 
!      coordinate direction. The Dirichlet condition is set for patches the 
!      normal to which is aligned with a coordinate direction.
!   2. This routine must be called after the boundary vertex normals are 
!      constructed, and after the boundary condition file is read (so that the
!      grid motion options for motion and smoothing are set).
!
! ******************************************************************************
!
! $Id: RFLU_SetMoveGridOptions.F90,v 1.15 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetMoveGridOptions(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError   
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch  
  USE ModMPI

  USE ModTools  
    
  IMPLICIT NONE

! **************************************I****************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================
 
  LOGICAL :: outsiderFlag,xFlatFlag,yFlatFlag,zFlatFlag
  CHARACTER(CHRLEN) :: moveBcTypeString,RCSIdentString  
  INTEGER :: errorFlag,ibv,iPatch
  INTEGER, DIMENSION(2) :: maxNormLoc
  REAL(RFREAL), PARAMETER :: eqTol = 1.0E-6_RFREAL
  REAL(RFREAL) :: term,xNorm,xNormMax,xNormMin,yNorm,yNormMax,yNormMin, &  
                  zNorm,zNormMax,zNormMin
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: globalVals,localVals  
  REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: bvnExt
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetMoveGridOptions.F90,v $ $Revision: 1.15 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetMoveGridOptions',&
  'RFLU_SetMoveGridOptions.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting grid motion options...'
    WRITE(STDOUT,'(A,3X,A,1X,E15.8)') SOLVER_NAME,'Equality tolerance:',eqTol
    WRITE(STDOUT,'(A,3X,A,3X,A,2X,A)') SOLVER_NAME,'Local','Global','Option'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate and initialize temporary memory
! ******************************************************************************

  ALLOCATE(bvnExt(MIN_VAL:MAX_VAL,XCOORD:ZCOORD,global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bvnExt')
  END IF ! global%error

  DO iPatch = 1,global%nPatches
    bvnExt(MIN_VAL,XCOORD:ZCOORD,iPatch) =  HUGE(1.0_RFREAL)
    bvnExt(MAX_VAL,XCOORD:ZCOORD,iPatch) = -HUGE(1.0_RFREAL)              
  END DO ! iPatch

  ALLOCATE(localVals(XCOORD:ZCOORD,global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals')
  END IF ! global%error 
  
  ALLOCATE(globalVals(XCOORD:ZCOORD,global%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals')
  END IF ! global%error 

! ******************************************************************************
! Compute local extrema of patch vertex-normal vectors
! ******************************************************************************
 
  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    xNormMin =  HUGE(1.0_RFREAL)
    xNormMax = -HUGE(1.0_RFREAL)
    yNormMin =  HUGE(1.0_RFREAL)
    yNormMax = -HUGE(1.0_RFREAL)
    zNormMin =  HUGE(1.0_RFREAL)
    zNormMax = -HUGE(1.0_RFREAL)    

    DO ibv = 1,pPatch%nBVertTot
      xNorm = pPatch%bvn(XCOORD,ibv) 
      yNorm = pPatch%bvn(YCOORD,ibv)
      zNorm = pPatch%bvn(ZCOORD,ibv)

      IF ( xNorm < xNormMin ) THEN 
        xNormMin = xNorm
      ELSE IF ( xNorm > xNormMax ) THEN 
        xNormMax = xNorm
      END IF ! xNorm

      IF ( yNorm < yNormMin ) THEN 
        yNormMin = yNorm
      ELSE IF ( yNorm > yNormMax ) THEN 
        yNormMax = yNorm
      END IF ! yNorm

      IF ( zNorm < zNormMin ) THEN 
        zNormMin = zNorm
      ELSE IF ( zNorm > zNormMax ) THEN 
        zNormMax = zNorm
      END IF ! zNorm
    END DO ! ibv

    bvnExt(MIN_VAL,XCOORD,pPatch%iPatchGlobal) = xNormMin
    bvnExt(MIN_VAL,YCOORD,pPatch%iPatchGlobal) = yNormMin
    bvnExt(MIN_VAL,ZCOORD,pPatch%iPatchGlobal) = zNormMin       

    bvnExt(MAX_VAL,XCOORD,pPatch%iPatchGlobal) = xNormMax
    bvnExt(MAX_VAL,YCOORD,pPatch%iPatchGlobal) = yNormMax
    bvnExt(MAX_VAL,ZCOORD,pPatch%iPatchGlobal) = zNormMax        
  END DO ! iPatch

! ******************************************************************************
! Compute global extrema of patch normal vectors
! ******************************************************************************

  DO iPatch = 1,global%nPatches
    localVals(XCOORD,iPatch) = bvnExt(MIN_VAL,XCOORD,iPatch)
    localVals(YCOORD,iPatch) = bvnExt(MIN_VAL,YCOORD,iPatch)
    localVals(ZCOORD,iPatch) = bvnExt(MIN_VAL,ZCOORD,iPatch) 
  END DO ! iPatch  

  CALL MPI_Reduce(localVals,globalVals,(ZCOORD-XCOORD+1)*global%nPatches, & 
                  MPI_RFREAL,MPI_MIN,MASTERPROC,global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  DO iPatch = 1,global%nPatches
    bvnExt(MIN_VAL,XCOORD,iPatch) = globalVals(XCOORD,iPatch)
    bvnExt(MIN_VAL,YCOORD,iPatch) = globalVals(YCOORD,iPatch)
    bvnExt(MIN_VAL,ZCOORD,iPatch) = globalVals(ZCOORD,iPatch)
        
    localVals(XCOORD,iPatch) = bvnExt(MAX_VAL,XCOORD,iPatch)
    localVals(YCOORD,iPatch) = bvnExt(MAX_VAL,YCOORD,iPatch)
    localVals(ZCOORD,iPatch) = bvnExt(MAX_VAL,ZCOORD,iPatch)                      
  END DO ! iPatch  
   
  CALL MPI_Reduce(localVals,globalVals,(ZCOORD-XCOORD+1)*global%nPatches, & 
                  MPI_RFREAL,MPI_MIN,MASTERPROC,global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  DO iPatch = 1,global%nPatches
    bvnExt(MAX_VAL,XCOORD,iPatch) = globalVals(XCOORD,iPatch)
    bvnExt(MAX_VAL,YCOORD,iPatch) = globalVals(YCOORD,iPatch)
    bvnExt(MAX_VAL,ZCOORD,iPatch) = globalVals(ZCOORD,iPatch)   
  END DO ! iPatch   

! ******************************************************************************
! Loop over patches to determine options on grid motion
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    pPatch%moveBcType = MOVEGRID_BCTYPE_NONE  

! ==============================================================================
!   Determine whether non-moving patches are flat
! ==============================================================================

    IF ( (pPatch%movePatch .EQV. .FALSE.) .AND. & 
         (pPatch%smoothGrid .EQV. .TRUE.) ) THEN 

! ------------------------------------------------------------------------------
!     Get extrema of normal vector
! ------------------------------------------------------------------------------
    
      xNormMin = bvnExt(MIN_VAL,XCOORD,pPatch%iPatchGlobal)
      xNormMax = bvnExt(MAX_VAL,XCOORD,pPatch%iPatchGlobal)

      yNormMin = bvnExt(MIN_VAL,YCOORD,pPatch%iPatchGlobal)        
      yNormMax = bvnExt(MAX_VAL,YCOORD,pPatch%iPatchGlobal)

      zNormMin = bvnExt(MIN_VAL,ZCOORD,pPatch%iPatchGlobal)
      zNormMax = bvnExt(MAX_VAL,ZCOORD,pPatch%iPatchGlobal)
    
! ------------------------------------------------------------------------------
!     Determine whether patch is flat by checking whether extrema have same 
!     sign (if yes, and extrema are equal, have flat patch; if no, and extrema 
!     are close to zero, also have flat patch)
! ------------------------------------------------------------------------------
                
! --- x-direction --------------------------------------------------------------    
        
      IF ( (SIGN(1.0_RFREAL,xNormMin) == SIGN(1.0_RFREAL,xNormMax)) ) THEN 
        IF ( FloatEqual(xNormMin,xNormMax,eqTol) .EQV. .TRUE. ) THEN 
          xFlatFlag = .TRUE.
        ELSE 
          xFlatFlag = .FALSE.
        END IF ! FloatEqual
      ELSE
         IF ( (FloatEqual(ABS(xNormMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
              (FloatEqual(ABS(xNormMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
          xFlatFlag = .TRUE.      
        ELSE
          xFlatFlag = .FALSE.
        END IF ! FloatEqual
      END IF ! FloatEqual

! --- y-direction --------------------------------------------------------------    

      IF ( (SIGN(1.0_RFREAL,yNormMin) == SIGN(1.0_RFREAL,yNormMax)) ) THEN 
        IF ( FloatEqual(yNormMin,yNormMax,eqTol) .EQV. .TRUE. ) THEN 
          yFlatFlag = .TRUE.
        ELSE 
          yFlatFlag = .FALSE.
        END IF ! FloatEqual
      ELSE
         IF ( (FloatEqual(ABS(yNormMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
              (FloatEqual(ABS(yNormMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
          yFlatFlag = .TRUE.      
        ELSE
          yFlatFlag = .FALSE.
        END IF ! FloatEqual
      END IF ! FloatEqual

! --- z-direction --------------------------------------------------------------    

      IF ( (SIGN(1.0_RFREAL,zNormMin) == SIGN(1.0_RFREAL,zNormMax)) ) THEN     
        IF ( FloatEqual(zNormMin,zNormMax,eqTol) .EQV. .TRUE. ) THEN 
          zFlatFlag = .TRUE.
        ELSE 
          zFlatFlag = .FALSE.
        END IF ! FloatEqual
      ELSE    
         IF ( (FloatEqual(ABS(zNormMin),EPSILON(1.0_RFREAL),eqTol)) .AND. & 
              (FloatEqual(ABS(zNormMax),EPSILON(1.0_RFREAL),eqTol)) ) THEN
          zFlatFlag = .TRUE.      
        ELSE
          zFlatFlag = .FALSE.
        END IF ! FloatEqual
      END IF ! FloatEqual
    
! ==============================================================================
!     Determine whether patch is flat and aligned with coordinate directions. 
!     If flat and aligned with coordinate directions, use Dirichlet conditions,
!     otherwise use Neumann conditions (but only if do not have any outsider
!     actual vertices).
! ==============================================================================
                  
      IF ( (xFlatFlag .EQV. .TRUE.) .AND. & 
           (yFlatFlag .EQV. .TRUE.) .AND. & 
           (zFlatFlag .EQV. .TRUE.) ) THEN 
        term = MAX(ABS(xNormMin),ABS(xNormMax), & 
                   ABS(yNormMin),ABS(yNormMax), &
                   ABS(zNormMin),ABS(zNormMax))             
                   
! ------------------------------------------------------------------------------
!       Patch is flat and aligned with coordinate axes
! ------------------------------------------------------------------------------
        
        IF ( FloatEqual(term,1.0_RFREAL) .EQV. .TRUE. ) THEN  
          IF ( term == ABS(xNormMin) .OR. term == ABS(xNormMax) ) THEN
            pPatch%moveBcType = XCOORD       
          ELSE IF ( term == ABS(yNormMin) .OR. term == ABS(yNormMax) ) THEN
            pPatch%moveBcType = YCOORD
          ELSE IF ( term == ABS(zNormMin) .OR. term == ABS(zNormMax) ) THEN
            pPatch%moveBcType = ZCOORD  
          END IF ! term                                  

! ------------------------------------------------------------------------------
!       Patch is flat but not aligned with coordinate axes, so set boundary 
!       condition to von Neumann
! ------------------------------------------------------------------------------        
        
        ELSE
          pPatch%moveBcType = MOVEGRID_BCTYPE_NEUMANN
        END IF ! xFlagFlag
      END IF ! xFlatFlag

! ==============================================================================
!     Check that patches have proper boundary condition. Patches with active 
!     smoothing must be flat, whether or not aligned with coordinate directions.
! ==============================================================================

      IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_NONE ) THEN 
        pPatch%smoothGrid = .FALSE.
         
        global%warnCounter = global%warnCounter + 1  
         
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_NONE ) THEN            
           WRITE(STDOUT,'(A,3X,A,I3,A)') SOLVER_NAME, &
                 '*** WARNING *** Invalid smoothing input for patch ',iPatch, &
                 '. Overriding user input for grid smoothing.'
        END IF ! global%myProcid                      
      END IF ! pPatch%smoothGrid
    END IF ! pPatch%movePatch

! ==============================================================================
!   Write out info
! ==============================================================================
  
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_NONE ) THEN 
        moveBcTypeString = 'None'   
      ELSE IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_NEUMANN ) THEN 
        moveBcTypeString = 'von Neumann (homogeneous)'         
      ELSE IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_DIRICHX ) THEN 
        moveBcTypeString = 'Dirichlet (homogeneous in x)'     
      ELSE IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_DIRICHY ) THEN 
        moveBcTypeString = 'Dirichlet (homogeneous in y)' 
      ELSE IF ( pPatch%moveBcType == MOVEGRID_BCTYPE_DIRICHZ ) THEN 
        moveBcTypeString = 'Dirichlet (homogeneous in z)' 
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! pPatch%moveBcType          
           
      WRITE(STDOUT,'(A,2X,I4,5X,I4,4X,A)') SOLVER_NAME,iPatch, & 
                                           pPatch%iPatchGlobal, & 
                                           TRIM(moveBcTypeString)
    END IF ! global%myProcid
  END DO ! iPatch

! ******************************************************************************
! Check settings - catch error before get to grid motion (caught there also)
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)
  
    IF ( (pPatch%smoothGrid .EQV. .TRUE.) .AND. & 
         (pPatch%moveBcType == MOVEGRID_BCTYPE_NONE) ) THEN 
      CALL ErrorStop(global,ERR_MOVEPATCH_BC_INVALID,__LINE__)
    END IF ! pPatch
  END DO ! iPatch

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(bvnExt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bvnExt')
  END IF ! global%error

  DEALLOCATE(localVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals')
  END IF ! global%error 
  
  DEALLOCATE(globalVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals')
  END IF ! global%error 

! ******************************************************************************
! End 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting grid motion options done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetMoveGridOptions

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetMoveGridOptions.F90,v $
! Revision 1.15  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.12  2005/04/15 15:07:25  haselbac
! Converted to MPI
!
! Revision 1.11  2004/10/19 19:29:24  haselbac
! Removed logic related to insider/outsider vertices
!
! Revision 1.10  2003/07/22 02:12:00  haselbac
! Added global%warnCounter
!
! Revision 1.9  2003/05/01 14:13:51  haselbac
! Added logic to deal with non-moving and non-smoothed patches
!
! Revision 1.8  2003/03/27 21:50:21  haselbac
! Fixed bug in determining outsiderFlag
!
! Revision 1.7  2003/03/15 18:57:22  haselbac
! Completed || of gm, added check for outsider vertices
!
! Revision 1.6  2003/02/20 20:15:23  haselbac
! Bug fix: nVertTot changed to nVert (from different working version)
!
! Revision 1.5  2003/02/20 19:49:20  haselbac
! Corrected bug in check
!
! Revision 1.4  2003/02/06 19:32:16  haselbac
! Bug fix (EPS instd of PREC), use tolerance, check for error
!
! Revision 1.3  2003/01/28 14:49:19  haselbac
! Only set bc for non-moving patches
!
! Revision 1.2  2002/11/27 20:27:17  haselbac
! Changed test for planarity and added output
!
! Revision 1.1  2002/11/26 15:28:56  haselbac
! Initial revision
!
! ******************************************************************************







