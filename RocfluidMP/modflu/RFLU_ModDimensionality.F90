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
! Purpose: Suite of routines related to two dimensional computations.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModDimensionality.F90,v 1.9 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModDimensionality

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_123D_CheckGeometryWrapper, & 
            RFLU_123D_CheckTopology
      
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModDimensionality.F90,v $ $Revision: 1.9 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  





! ******************************************************************************
!
! Purpose: Check whether geometry is correct for 1d and 2d runs.
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

  SUBROUTINE RFLU_123D_CheckGeometryWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    TYPE(t_global), POINTER :: global      
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_123D_CheckGeometryWrapper',&
  'RFLU_ModDimensionality.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Check geometry for 2d runs
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        CALL RFLU_123D_CheckGeometryKernel(pRegion,YCOORD)
        CALL RFLU_123D_CheckGeometryKernel(pRegion,ZCOORD)      
      CASE ( 2 ) 
        CALL RFLU_123D_CheckGeometryKernel(pRegion,ZCOORD)
      CASE ( 3 ) ! Defensive coding
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_123D_CheckGeometryWrapper








! ******************************************************************************
!
! Purpose: Check whether geometry is correct for 1d and 2d runs.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   dir         Coordinate direction
!
! Output: None.
!
! Notes: 
!   1. dir-component of normals MUST be close to machine precision for all 
!      patches except virtual ones, for which they must be close to unity.
!
! ******************************************************************************

  SUBROUTINE RFLU_123D_CheckGeometryKernel(pRegion,dir)

    USE ModTools, ONLY: FloatEqual

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    INTEGER, INTENT(IN) :: dir
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: iPatch
    REAL(RFREAL) :: ndMax,nxMax,ndMin,nxMin,nTol
    TYPE(t_grid), POINTER :: pGrid      
    TYPE(t_patch), POINTER :: pPatch            
    TYPE(t_global), POINTER :: global      

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_123D_CheckGeometryKernel',&
  'RFLU_ModDimensionality.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking geometry...'
      WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Component:',dir
    END IF ! global%myProcid

    pGrid => pRegion%grid

    nTol =  1.0E-14_RFREAL

! ******************************************************************************
!   Interior faces: dir-component of normal must be close to zero
! ******************************************************************************

    IF ( pGrid%nFacesTot > 0 ) THEN 
      ndMin = MINVAL(pGrid%fn(dir,1:pGrid%nFacesTot))
      ndMax = MAXVAL(pGrid%fn(dir,1:pGrid%nFacesTot))      

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH) THEN     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Extrema of face-normal vectors:'
        WRITE(STDOUT,'(A,5X,A,1X,E23.16)') SOLVER_NAME,'Tolerance:',nTol
        WRITE(STDOUT,'(A,7X,A,1X,2(1X,E23.16))') SOLVER_NAME,'Interior:', &
                                                 ndMin,ndMax
      END IF ! global%myProcid

      IF ( (ABS(ndMin) > nTol) .OR. (ABS(ndMax) > nTol) ) THEN 
        CALL ErrorStop(global,ERR_FACE_NORMAL_INVALID,__LINE__)
      END IF ! ABS(ndMin)
  
      IF ( (FloatEqual(ABS(ndMin),0.0_RFREAL,nTol) .EQV. .FALSE.) .OR. & 
           (FloatEqual(ABS(ndMax),0.0_RFREAL,nTol) .EQV. .FALSE.) ) THEN 
        CALL ErrorStop(global,ERR_FACE_NORMAL_INVALID,__LINE__)
      END IF ! FloatEqual     
    END IF ! pGrid%nFacesTot

! ******************************************************************************
!   Boundary faces: dir-component of normal must be close to zero for faces on
!   non-virtual patches, otherwise must be close to unity.
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%nBFacesTot > 0 ) THEN 
        ndMin = MINVAL(pPatch%fn(dir,1:pPatch%nBFacesTot))
        ndMax = MAXVAL(pPatch%fn(dir,1:pPatch%nBFacesTot))  

        nxMin = MINVAL(pPatch%fn(XCOORD,1:pPatch%nBFacesTot))
        nxMax = MAXVAL(pPatch%fn(XCOORD,1:pPatch%nBFacesTot))  

        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel >= VERBOSE_HIGH) THEN
          WRITE(STDOUT,'(A,7X,A,1X,I2,A,1X,2(1X,E23.16))') & 
                SOLVER_NAME,'Patch',iPatch,':',ndMin,ndMax
        END IF ! global%myProcid

! ==============================================================================
!       Non-virtual patches: dir-component of normal must be close to zero 
! ==============================================================================

        IF ( pPatch%bcType /= BC_VIRTUAL ) THEN        
          IF ( (FloatEqual(ABS(ndMin),0.0_RFREAL,nTol) .EQV. .FALSE.) .OR. & 
               (FloatEqual(ABS(ndMax),0.0_RFREAL,nTol) .EQV. .FALSE.) ) THEN 
            CALL ErrorStop(global,ERR_FACE_NORMAL_INVALID,__LINE__)
          END IF ! FloatEqual 
  
! ==============================================================================
!       Virtual patches: x-component of normal must be close to zero
! ==============================================================================

        ELSE 
          IF ( (FloatEqual(ABS(nxMin),0.0_RFREAL,nTol) .EQV. .FALSE.) .OR. & 
               (FloatEqual(ABS(nxMax),0.0_RFREAL,nTol) .EQV. .FALSE.) ) THEN 
            CALL ErrorStop(global,ERR_FACE_NORMAL_INVALID,__LINE__)
          END IF ! FloatEqual     
        END IF ! pPatch%bcType            
      END IF ! pPatch%nBFacesTot
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_123D_CheckGeometryKernel






! ******************************************************************************
!
! Purpose: Check whether topology is correct.
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

  SUBROUTINE RFLU_123D_CheckTopology(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: icl,ifg,ifl,ifl2,iPatch,iPatchCntr
    TYPE(t_grid), POINTER :: pGrid      
    TYPE(t_patch), POINTER :: pPatch            
    TYPE(t_global), POINTER :: global      

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_123D_CheckTopology',&
  'RFLU_ModDimensionality.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking topology...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Check dimensionality
! ******************************************************************************

! ==============================================================================
!   Check that have only appropriate cells     
! ==============================================================================
    
    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        IF ( pGrid%nTetsTot /= 0 .OR. &
             pGrid%nPrisTot /= 0 .OR. &
             pGrid%nPyrsTot /= 0 ) THEN
          CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
        END IF ! pGrid%nTetsTot    
      CASE ( 2 ) 
        IF ( pGrid%nTetsTot /= 0 .OR. pGrid%nPyrsTot /= 0 ) THEN
          CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
        END IF ! pGrid%nTetsTot
      CASE ( 3 ) 
      CASE DEFAULT ! Defensive coding
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens
    
! ==============================================================================
!   Must have specified number of virtual patches    
! ==============================================================================

    iPatchCntr = 0 

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch) 

      IF ( pPatch%bcType == BC_VIRTUAL ) THEN 
        iPatchCntr = iPatchCntr + 1
      END IF ! pPatch%bcType
    END DO ! iPatch

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 1 )         
        IF ( iPatchCntr /= 4 ) THEN 
          CALL ErrorStop(global,ERR_NUM_BC_VIRTUAL,__LINE__)        
        END IF ! iPatchCntr    
      CASE ( 2 )         
        IF ( iPatchCntr /= 2 ) THEN 
          CALL ErrorStop(global,ERR_NUM_BC_VIRTUAL,__LINE__)        
        END IF ! iPatchCntr
      CASE ( 3 )
        IF ( iPatchCntr /= 0 ) THEN 
          CALL ErrorStop(global,ERR_NUM_BC_VIRTUAL,__LINE__)        
        END IF ! iPatchCntr       
      CASE DEFAULT ! Defensive coding
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens    

! ==============================================================================
!   Each cell must have specified number of boundary faces. NOTE require 
!   cell-to-face connectivity array. For prisms, boundary faces must be 
!   triangular faces. 
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%dimens )     

! ------------------------------------------------------------------------------
!     One dimension
! ------------------------------------------------------------------------------      

      CASE ( 1 ) 
      
! ----- Hexahedra --------------------------------------------------------------
      
        DO icl = 1,pGrid%nHexsTot
          iPatchCntr = 0 
        
          DO ifl = 1,6
            iPatch = pGrid%hex2f(1,ifl,icl)
            
            IF ( iPatch > 0 ) THEN
              pPatch => pRegion%patches(iPatch)
             
              IF ( pPatch%bcType == BC_VIRTUAL ) THEN            
                iPatchCntr = iPatchCntr + 1
              END IF ! pPatch%bcType
            END IF ! iPatch
          END DO ! ifl
          
          IF ( iPatchCntr /= 4 ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! iPatchCntr
        END DO ! icl

! ------------------------------------------------------------------------------
!     Two dimensions
! ------------------------------------------------------------------------------      

      CASE ( 2 ) 
      
! ----- Hexahedra --------------------------------------------------------------
      
        DO icl = 1,pGrid%nHexsTot
          iPatchCntr = 0 
        
          DO ifl = 1,6
            iPatch = pGrid%hex2f(1,ifl,icl)
            
            IF ( iPatch > 0 ) THEN
              pPatch => pRegion%patches(iPatch)
             
              IF ( pPatch%bcType == BC_VIRTUAL ) THEN            
                iPatchCntr = iPatchCntr + 1
              END IF ! pPatch%bcType
            END IF ! iPatch
          END DO ! ifl
          
          IF ( iPatchCntr /= 2 ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! iPatchCntr
        END DO ! icl

! ----- Prisms -----------------------------------------------------------------
        
        DO icl = 1,pGrid%nPrisTot
          iPatchCntr = 0 
        
          DO ifl = 1,5
            iPatch = pGrid%pri2f(1,ifl,icl)
            ifl2   = pGrid%pri2f(2,ifl,icl)
            
            IF ( iPatch > 0 ) THEN
              pPatch => pRegion%patches(iPatch)
             
              IF ( pPatch%bcType == BC_VIRTUAL ) THEN            
                iPatchCntr = iPatchCntr + 1
                                
                IF ( pPatch%bf2v(4,ifl2) /= VERT_NONE ) THEN 
                  CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
                END IF ! pGrid%f2v                
              END IF ! pPatch%bcType
            END IF ! iPatch            
          END DO ! ifl
          
          IF ( iPatchCntr /= 2 ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! iPatchCntr
        END DO ! icl       

! ------------------------------------------------------------------------------
!     Three dimensions and default
! ------------------------------------------------------------------------------      

      CASE ( 3 ) 
      CASE DEFAULT ! Defensive coding
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking topology done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_123D_CheckTopology




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModDimensionality


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDimensionality.F90,v $
! Revision 1.9  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2007/03/07 03:20:05  haselbac
! Added IFs on nFacesTot to avoid problems with single-cell grid
!
! Revision 1.6  2007/02/27 13:03:55  haselbac
! Enabled 1d computations
!
! Revision 1.5  2005/12/11 15:55:03  haselbac
! Bug fix: Added missing IF on myProcid
!
! Revision 1.4  2005/11/09 01:23:05  haselbac
! Increased tolerance again, did not increase version number
!
! Revision 1.3  2005/11/09 01:18:20  haselbac
! Changed tolerance after finding problems with Manojs cylinder comp
!
! Revision 1.2  2005/11/04 14:06:48  haselbac
! Renamed existing routine, added new routine to check geom
! 
! ******************************************************************************
  









