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
! Purpose: Suite of routines to carry out split tree operations.
!
! Description: None.
!
! Notes: To create and use a split tree, one has to take the following steps:
!   1. 
!   2. 
!   3. Deallocate the table by calling DestroyHashTable.
!
!*******************************************************************************
!
! $Id: RFLU_ModSplitTree.F90,v 1.6 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModSplitTree

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_CreateSplitTree, & 
            RFLU_BuildSplitTree, & 
            RFLU_DestroySplitTree

  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
  INTEGER, PARAMETER :: NLEVELS_MAX  = 50,    & ! Must be greater than 1
                        NBUCKETS_MAX = 10000, & ! Must be greater than 2
                        NPOINTS_MAX  = 5       
  
  INTEGER, PRIVATE :: nPoints
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: pointList     
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: tree     
       
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: pointXyz     
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
! ==============================================================================
!   Create split tree
! ==============================================================================  

    SUBROUTINE RFLU_CreateSplitTree(global,nDataPoints)
    
      IMPLICIT NONE   
        
      INTEGER, INTENT(IN) :: nDataPoints 
      TYPE(t_global), POINTER :: global  
        
      INTEGER :: errorFlag,ip         
        
      CALL RegisterFunction(global,'RFLU_CreateSplitTree',&
  'RFLU_ModSplitTree.F90')

! ------------------------------------------------------------------------------
!     Copy argument into nPoints variable
! ------------------------------------------------------------------------------ 

      nPoints = nDataPoints

! ------------------------------------------------------------------------------
!     Allocate memory
! ------------------------------------------------------------------------------ 

      ALLOCATE(tree(7,NBUCKETS_MAX),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error)

      ALLOCATE(pointList(nPoints),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

      ALLOCATE(pointXyz(nPoints),STAT=errorFlag)
      global%error = errorFlag  
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

! ------------------------------------------------------------------------------
!     Initialize
! ------------------------------------------------------------------------------ 

      tree(:,:) = 0

      DO ip = 1,nPoints
        pointList(ip) = ip
      END DO ! ip
    
      CALL DeregisterFunction(global)    

    END SUBROUTINE RFLU_CreateSplitTree

! ==============================================================================
!   Build split tree
! ==============================================================================  

    SUBROUTINE RFLU_BuildSplitTree(global,xyz)

      IMPLICIT NONE
      
      INTEGER :: errorFlag,iBranch1,iBranch2,iBucket,iBucketLast,ipl,ipm,ipg, & 
                 iSplitDir,nBuckets,nLevels
      
      REAL(RFREAL), INTENT(IN) :: xyz(3,nPoints)
      TYPE(t_global), POINTER :: global
      
! ------------------------------------------------------------------------------
!     Start
! ------------------------------------------------------------------------------      
      
      CALL RegisterFunction(global,'RFLU_BuildSplitTree',&
  'RFLU_ModSplitTree.F90')      

      nLevels  = 1
      nBuckets = 1
      
      iBucket     = 1
      iBucketLast = 1
      iSplitDir   = 1
      
      tree(1,iBucket) = 1
      tree(2,iBucket) = nPoints
      tree(3,iBucket) = nPoints
      tree(4,iBucket) = iSplitDir
      
! ------------------------------------------------------------------------------
!     Loop over levels
! ------------------------------------------------------------------------------      
      
      OUTER: DO
        iSplitDir = iSplitDir + 1
        IF ( iSplitDir > 3 ) THEN 
          iSplitDir = 1
        END IF ! iSplitDir
      
! ----- Loop over buckets      
      
        INNER: DO 
     
! ------- Determine whether bucket should be split     
           
          IF ( tree(3,iBucket) > NPOINTS_MAX ) THEN ! Split bucket

! --------- Check whether more buckets can be generated

            IF ( nBuckets > NBUCKETS_MAX-2 ) THEN 
              EXIT OUTER
            END IF ! nBuckets
          
! --------- Determine where bucket should be split - use median 
          
            ipl = 0
            DO ipg = tree(1,iBucket),tree(2,iBucket)
              ipl = ipl + 1
              pointXyz(ipl) = xyz(iSplitDir,pointList(ipg))
            END DO ! ipg
            
            CALL QuickSortRFREALInteger(pointXyz(1:tree(3,iBucket)), &
                          pointList(tree(1,iBucket):tree(2,iBucket)), & 
                          tree(3,iBucket))
          
            ipm = tree(3,iBucket)/2 + 1 ! NOTE integer division
          
! --------- Update information for current bucket          
          
            tree(5,iBucket) = pointList(ipm)          
          
            iBranch1 = iBucket + 1
            iBranch2 = iBucket + 2
            
            tree(6,iBucket) = iBranch1
            tree(7,iBucket) = iBranch2          
          
! --------- Create new buckets

            nBuckets = nBuckets + 2
            
            tree(1,iBranch1) = tree(1,iBucket)
            tree(2,iBranch1) = ipm - 1
            tree(3,iBranch1) = tree(2,iBranch1) - tree(1,iBranch1) + 1        
            tree(4,iBranch1) = iSplitDir
          
            tree(1,iBranch2) = ipm
            tree(2,iBranch2) = tree(3,iBucket)
            tree(3,iBranch2) = tree(2,iBranch2) - tree(1,iBranch2) + 1
            tree(4,iBranch2) = tree(4,iBranch1)           
          END IF ! tree(3,iBucket)  
          
! ------- Check whether more buckets to be split          
        
          IF ( iBucket < nBuckets ) THEN 
            IF ( iBucket /= iBucketLast ) THEN
              iBucket = iBucket + 1             
            ELSE 
              iBucketLast = nBuckets
              EXIT INNER
            END IF ! iBucket
          ELSE
            EXIT OUTER  
          END IF ! iBucket
        END DO INNER 
      
! ----- Check whether more levels can be inserted      
      
        IF ( nLevels < NLEVELS_MAX ) THEN 
          nLevels   = nLevels   + 1
          iBucket   = iBucket   + 1
        ELSE  
          EXIT OUTER
        END IF ! nLevels
      END DO OUTER

      WRITE(*,*) nLevels,nBuckets

      DO ipl = 1,nBuckets
        WRITE(*,*) ipl,tree(1:7,ipl)
      END DO ! ipl

      CALL DeregisterFunction(global)  

    END SUBROUTINE RFLU_BuildSplitTree

! ==============================================================================
!   Destroy split tree
! ==============================================================================  
   
    SUBROUTINE RFLU_DestroySplitTree(global)
    
      IMPLICIT NONE

      TYPE(t_global), POINTER :: global

      INTEGER :: errorFlag

      CALL RegisterFunction(global,'RFLU_DestroySplitTree',&
  'RFLU_ModSplitTree.F90')

      DEALLOCATE(tree,STAT=errorFlag)
      global%error = errorFlag  
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error)

      DEALLOCATE(pointList,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error

      DEALLOCATE(pointXyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error

      CALL DeregisterFunction(global)

    END SUBROUTINE RFLU_DestroySplitTree


  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModSplitTree


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModSplitTree.F90,v $
!   Revision 1.6  2008/12/06 08:44:24  mtcampbe
!   Updated license.
!
!   Revision 1.5  2008/11/19 22:17:35  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.4  2004/01/22 16:03:59  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.3  2002/10/08 15:49:21  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.2  2002/09/09 15:12:12  haselbac
!   global now under regions
!
!   Revision 1.1  2002/04/11 18:48:48  haselbac
!   Initial revision
!
!
! ******************************************************************************









