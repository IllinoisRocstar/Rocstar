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
! Purpose: Pick special faces.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region data
! 
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PickSpecialFaces.F90,v 1.3 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PickSpecialFaces(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModSortSearch

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

  CHARACTER :: infoType,stencilType
  CHARACTER(CHRLEN) :: RCSIdentString  
  INTEGER :: errorFlag,faceIndx,iFacesSpecial,patchIndx
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickSpecialFaces.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PickSpecialFaces', &
                        'RFLU_PickSpecialFaces.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special faces...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal  
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and initialize
! ******************************************************************************

  pGrid => pRegion%grid

  iFacesSpecial = 0 
  pGrid%facesSpecial(1:2,1:NFACES_SPECIAL_MAX) = 0

! ******************************************************************************
! Get information from user
! ******************************************************************************

  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter information on special faces:' 
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'b - boundary face'    
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'i - interior face'      
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'q - quit'

! ******************************************************************************
! Set up infinite loop
! ******************************************************************************

  DO
  
! ==============================================================================
!   Enter information type
! ==============================================================================
   
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
    READ(STDIN,'(A)') infoType
  
    SELECT CASE ( infoType )

! ------------------------------------------------------------------------------    
!     Boundary face
! ------------------------------------------------------------------------------

      CASE ( 'b' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter patch index:'
        READ(STDIN,*,IOSTAT=errorFlag) patchIndx

        IF ( errorFlag /= ERR_NONE ) THEN 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag        
          
        IF ( patchIndx > 0 .AND. patchIndx <= pGrid%nPatches ) THEN
          pPatch => pRegion%patches(patchIndx)
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter face index:'
          READ(STDIN,*,IOSTAT=errorFlag) faceIndx           

          IF ( errorFlag /= ERR_NONE ) THEN
            global%warnCounter = global%warnCounter + 1           
           
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.'
            CYCLE
          END IF ! errorFlag 

          IF ( faceIndx > 0 .AND. faceIndx <= pPatch%nBFacesTot ) THEN
            IF ( iFacesSpecial == NFACES_SPECIAL_MAX ) THEN 
              CALL ErrorStop(global,ERR_NFACES_SPECIAL_MAX,__LINE__)
            END IF ! iFacesSpecial         
           
            iFacesSpecial = iFacesSpecial + 1          
            pGrid%facesSpecial(1,iFacesSpecial) = patchIndx
            pGrid%facesSpecial(2,iFacesSpecial) = faceIndx            
  
            WRITE(STDOUT,'(A,5X,A,1X,I8)') SOLVER_NAME,'Added face:',faceIndx        
          ELSE 
            global%warnCounter = global%warnCounter + 1           
          
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.' 
            CYCLE 
          END IF ! faceIndx        
        ELSE 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! patchIndx
   
! ------------------------------------------------------------------------------
!     Interior face     
! ------------------------------------------------------------------------------

      CASE ( 'i' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter interior face index:'
        READ(STDIN,*,IOSTAT=errorFlag) faceIndx
 
        IF ( errorFlag /= ERR_NONE ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag         
        
        IF ( faceIndx > 0 .AND. faceIndx <= pGrid%nFacesTot ) THEN
          IF ( iFacesSpecial == NFACES_SPECIAL_MAX ) THEN ! NOTE
            CALL ErrorStop(global,ERR_NFACES_SPECIAL_MAX,__LINE__)
          END IF ! iFacesSpecial       
         
          iFacesSpecial = iFacesSpecial + 1          
          pGrid%facesSpecial(1,iFacesSpecial) = 0
          pGrid%facesSpecial(2,iFacesSpecial) = faceIndx
                                        
          WRITE(STDOUT,'(A,5X,A,1X,I8)') SOLVER_NAME,'Added face:',faceIndx   
        ELSE 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! faceIndx        
            
! ------------------------------------------------------------------------------
!     Quit       
! ------------------------------------------------------------------------------
      
      CASE ( 'q' )
        EXIT
        
! ------------------------------------------------------------------------------
!     Default        
! ------------------------------------------------------------------------------
        
      CASE DEFAULT 
        global%warnCounter = global%warnCounter + 1 
      
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
        CYCLE 
    END SELECT  
  END DO ! <empty> 

! ******************************************************************************
! Set number of special faces
! ******************************************************************************

  pGrid%nFacesSpecial = iFacesSpecial

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special faces done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickSpecialFaces

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickSpecialFaces.F90,v $
! Revision 1.3  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/27 02:02:23  haselbac
! Initial revision
!
! ******************************************************************************







