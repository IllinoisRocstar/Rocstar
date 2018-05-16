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
! Purpose: Suite of routines to build boundary-face gradient access list.
!
! Description: None.
!
! Notes: 
!   1. The routine which builds the access lists MUST be called after the face 
!      lists have been built. 
!
!*******************************************************************************
!
! $Id: RFLU_ModBFaceGradAccessList.F90,v 1.7 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModBFaceGradAccessList

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_NullifyBFaceGradAccessList, & 
            RFLU_CreateBFaceGradAccessList, & 
            RFLU_BuildBFaceGradAccessList, & 
            RFLU_DestroyBFaceGradAccessList
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModBFaceGradAccessList.F90,v $ $Revision: 1.7 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS
  


! ******************************************************************************
!   Nullify boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_NullifyBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_NullifyBFaceGradAccessList',&
  'RFLU_ModBFaceGradAccessList.F90') 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Nullifying boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
          
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

! TEMPORARY : removing usage of bf2bg from everywhere
!        NULLIFY(pPatch%bf2bg)
!        NULLIFY(pPatch%bf2bgTot)
      END DO ! iPatch          
        
! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Nullifying boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_NullifyBFaceGradAccessList




! ******************************************************************************
!   Create boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_CreateBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_CreateBFaceGradAccessList',&
  'RFLU_ModBFaceGradAccessList.F90') 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Creating boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
          
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

! TEMPORARY : removing usage of bf2bg from everywhere
!        ALLOCATE(pPatch%bf2bg(BF2BG_BEG:BF2BG_END),STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2bg')          
!        END IF ! global%error
        
! TEMPORARY : removing usage of bf2bg from everywhere
!        ALLOCATE(pPatch%bf2bgTot(BF2BG_BEG:BF2BG_END),STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2bgTot')          
!        END IF ! global%error        
      END DO ! iPatch        
        
! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Creating boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_CreateBFaceGradAccessList
    



  
! ******************************************************************************
!   Build boundary-face gradient access list
! ****************************************************************************** 
    
    SUBROUTINE RFLU_BuildBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================

      TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!     Locals
! ==============================================================================

      INTEGER :: bf2bgLast,bf2bgTotLast,iPatch
      TYPE(t_global), POINTER :: global
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!     Start
! ******************************************************************************

      global => pRegion%global

      CALL RegisterFunction(global,'RFLU_BuildBFaceGradAccessList',&
  'RFLU_ModBFaceGradAccessList.F90')

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN    
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Building boundary-face gradient access lists...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid  

! ==============================================================================
!     Loop over patches
! ==============================================================================

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( iPatch /= 1 ) THEN    
! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bg(BF2BG_BEG) = bf2bgLast + 1        
!          pPatch%bf2bg(BF2BG_END) = pPatch%bf2bg(BF2BG_BEG) & 
!                                  + pPatch%nBFaces - 1

! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bgTot(BF2BG_BEG) = bf2bgTotLast + 1        
!          pPatch%bf2bgTot(BF2BG_END) = pPatch%bf2bgTot(BF2BG_BEG) & 
!                                     + pPatch%nBFacesTot - 1
        ELSE             
! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bg(BF2BG_BEG) = 1 
!          pPatch%bf2bg(BF2BG_END) = pPatch%nBFaces
 
! TEMPORARY : removing usage of bf2bg from everywhere         
!          pPatch%bf2bgTot(BF2BG_BEG) = 1 
!          pPatch%bf2bgTot(BF2BG_END) = pPatch%nBFacesTot              
        END IF ! iPatch

! TEMPORARY : removing usage of bf2bg from everywhere
!        bf2bgLast    = pPatch%bf2bg(BF2BG_END) 
!        bf2bgTotLast = pPatch%bf2bgTot(BF2BG_END)            
      END DO ! iPatch

! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Building boundary-face gradient access lists done.'
      END IF ! global%verbLevel

      CALL DeregisterFunction(global)

    END SUBROUTINE RFLU_BuildBFaceGradAccessList
    
      
      
  
 


! ******************************************************************************
!   Destroy boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_DestroyBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_DestroyBFaceGradAccessList',&
  'RFLU_ModBFaceGradAccessList.F90') 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Destroying boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
      
! TEMPORARY : removing usage of bf2bg from everywhere
!        DEALLOCATE(pPatch%bf2bg,STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2bg')          
!        END IF ! global%error
 
! TEMPORARY : removing usage of bf2bg from everywhere       
!        DEALLOCATE(pPatch%bf2bgTot,STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2bgTot')          
!        END IF ! global%error        
      END DO ! iPatch          
        
! ******************************************************************************
!     Nullify memory
! ******************************************************************************

      CALL RFLU_NullifyBFaceGradAccessList(pRegion)

! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Destroying boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_DestroyBFaceGradAccessList
 




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModBFaceGradAccessList


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModBFaceGradAccessList.F90,v $
!   Revision 1.7  2008/12/06 08:44:20  mtcampbe
!   Updated license.
!
!   Revision 1.6  2008/11/19 22:17:31  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.5  2006/08/19 15:39:21  mparmar
!   Removed bf2bg,bf2bgTot
!
!   Revision 1.4  2004/05/25 01:34:16  haselbac
!   Added code for bf2bgTot array; needed for Rocturb Sij access
!
!   Revision 1.3  2004/01/22 16:03:58  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.2  2003/12/09 03:58:01  haselbac
!   Bug fix for building of bf2bg list (no virtual faces)
!
!   Revision 1.1  2003/11/03 03:51:54  haselbac
!   Initial revision
!
! ******************************************************************************










