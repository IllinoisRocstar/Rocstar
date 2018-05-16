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
! Purpose: Collection of routines to write grids in STL format.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModSTL.F90,v 1.4 2008/12/06 08:44:55 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModSTL

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModSTL.F90,v $ $Revision: 1.4 $'        

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_STL_WriteSurfGridASCII

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  



! ******************************************************************************
!
! Purpose: Write surface grid in ASCII STL format.
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

  SUBROUTINE RFLU_STL_WriteSurfGridASCII(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: writePatchFlag
    CHARACTER(1) :: writePatchChar 
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,ifl,iFile,iPatch,ivg,ivl,ivl2
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_STL_WriteSurfGridASCII', &
                          'RFLU_ModSTL.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing surface grid in ASCII STL format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set pointers and variables
! ==============================================================================

    pGrid => pRegion%grid

! ******************************************************************************
!   Writing file
! ******************************************************************************

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.stl',iFileName)  

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Write file by looping over patches
! ==============================================================================

    WRITE(iFile,'(A)') 'solid'

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!     Determine whether patch is to be written
! ------------------------------------------------------------------------------      
      
      WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME,'Write patch', & 
                                          iPatch,'to STL file? (y/n)'
      READ(STDIN,'(A)',IOSTAT=errorFlag) writePatchChar 
      
      IF ( errorFlag == ERR_NONE ) THEN 
        IF ( writePatchChar == 'y' ) THEN 
	  writePatchFlag = .TRUE.
	ELSE
	  writePatchFlag = .FALSE. 
	END IF ! choice
      ELSE 
        writePatchFlag = .FALSE.
      END IF ! errorFlag

! ------------------------------------------------------------------------------
!     Write patch. NOTE quadrilateral faces written as two triangles.
! ------------------------------------------------------------------------------      
      
      IF ( writePatchFlag .EQV. .TRUE. ) THEN 
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A)') SOLVER_NAME,'Writing patch', & 
                                              iPatch,'to STL file...'      
        END IF ! global%verbLevel 
      
! ----- Triangles --------------------------------------------------------------
      
	DO ifl = 1,pPatch%nBTris
          WRITE(iFile,'(A,3(1X,E13.6))') 'facet normal', & 
	                        	 pPatch%fn(XCOORD:ZCOORD,ifl)                                       
          WRITE(iFile,'(A)') 'outer loop'                                         

          DO ivl = 1,3		
            ivg = pPatch%bTri2v(ivl,ifl)

            WRITE(iFile,'(A,3(1X,E13.6))') 'vertex',pGrid%xyz(XCOORD:ZCOORD,ivg)
          END DO ! ivl   

          WRITE(iFile,'(A)') 'endloop'          
          WRITE(iFile,'(A)') 'endfacet'                                              				       
	END DO ! ifl  

! ----  Quadrilaterals ---------------------------------------------------------
      
	DO ifl = 1,pPatch%nBQuads

! ------- Triangle 1       

          WRITE(iFile,'(A,3(1X,E13.6))') 'facet normal', & 
	                        	 pPatch%fn(XCOORD:ZCOORD,ifl)                                       
          WRITE(iFile,'(A)') 'outer loop'                                         

          DO ivl = 1,3		
            ivg = pPatch%bQuad2v(ivl,ifl)

            WRITE(iFile,'(A,3(1X,E13.6))') 'vertex',pGrid%xyz(XCOORD:ZCOORD,ivg)
          END DO ! ivl   

          WRITE(iFile,'(A)') 'endloop'          
          WRITE(iFile,'(A)') 'endfacet'                  

! ------- Triangle 2       

          WRITE(iFile,'(A,3(1X,E13.6))') 'facet normal', & 
	                        	 pPatch%fn(XCOORD:ZCOORD,ifl)                                       
          WRITE(iFile,'(A)') 'outer loop'         

          DO ivl = 3,5
            IF ( ivl == 5 ) THEN 
              ivl2 = 1
            ELSE 
              ivl2 = ivl
            END IF ! ivl 

            ivg = pPatch%bQuad2v(ivl2,ifl)

            WRITE(iFile,'(A,3(1X,E13.6))') 'vertex',pGrid%xyz(XCOORD:ZCOORD,ivg)
          END DO ! ivl         

          WRITE(iFile,'(A)') 'endloop'          
          WRITE(iFile,'(A)') 'endfacet'                                              				       
	END DO ! ifl
	
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A)') SOLVER_NAME,'Writing patch', & 
                                              iPatch,'to STL file done.'      
        END IF ! global%verbLevel 	
      END IF ! writePatchFlag          
    END DO ! iPatch

    WRITE(iFile,'(A)') 'endsolid'

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing surface grid in ASCII STL format done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_STL_WriteSurfGridASCII
  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModSTL

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModSTL.F90,v $
! Revision 1.4  2008/12/06 08:44:55  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/07/07 18:01:14  haselbac
! Bug fix (found on alc) for loop index
!
! Revision 1.1  2005/07/07 03:51:55  haselbac
! Initial revision
!
! ******************************************************************************
  
  







