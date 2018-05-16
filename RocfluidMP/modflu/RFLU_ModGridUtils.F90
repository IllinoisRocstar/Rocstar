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
! Purpose: Collection of utility routines for grid.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGridUtils.F90,v 1.3 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGridUtils

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGridUtils.F90,v $ $Revision: 1.3 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_DistortGrid

! ==============================================================================
! Private functions
! ==============================================================================



! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
     
   
   
! ******************************************************************************
!
! Purpose: Distort grid.
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
! ******************************************************************************

SUBROUTINE RFLU_DistortGrid(pRegion)
  
  USE ModRandom, ONLY: Rand1Uniform
  
  USE RFLU_ModBoundaryTests, ONLY: RFLU_TestIsBoundaryVertex    
      
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

  INTEGER :: errorFlag,ivg
  REAL(RFREAL) :: dx,dy,dz,xRand,yRand,zRand
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DistortGrid',&
  'RFLU_ModGridUtils.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN       
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Distorting grid...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  dx = global%distortX
  dy = global%distortY
  dz = global%distortZ

! ******************************************************************************
! Loop over vertices
! ******************************************************************************

  DO ivg = 1,pGrid%nVertTot
    IF ( RFLU_TestIsBoundaryVertex(pRegion,ivg) .EQV. .FALSE. ) THEN 
      xRand = Rand1Uniform(pRegion%randData)
      yRand = Rand1Uniform(pRegion%randData)
      zRand = Rand1Uniform(pRegion%randData)            
    
      pGrid%xyz(XCOORD,ivg) = pGrid%xyz(XCOORD,ivg) + xRand*dx
      pGrid%xyz(YCOORD,ivg) = pGrid%xyz(YCOORD,ivg) + yRand*dy
      pGrid%xyz(ZCOORD,ivg) = pGrid%xyz(ZCOORD,ivg) + zRand*dz           
    END IF ! RFLU_TestIsBoundaryVertex
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN       
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Distorting grid done ..'
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DistortGrid
   
   
   
  


END MODULE RFLU_ModGridUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGridUtils.F90,v $
! Revision 1.3  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/04 02:59:30  haselbac
! Initial revision
!
! ******************************************************************************









