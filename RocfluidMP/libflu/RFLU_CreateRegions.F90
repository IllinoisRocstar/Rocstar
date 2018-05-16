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
!******************************************************************************
!
! Purpose: Allocate memory for region pointer.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!   iLev                Level index
!   levels              Pointer to levels
!
! Output: None. 
!
! Notes: 
!   1. This routine exists only for speeding up compilation with optimization.
!      See routine RFLU_BuildDataStruct for detailed comments.
!
!******************************************************************************
!
! $Id: RFLU_CreateRegions.F90,v 1.3 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_CreateRegions(global,iLev,levels)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModError

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iLev
  TYPE(t_global), POINTER :: global
  TYPE(t_level), DIMENSION(:), POINTER :: levels

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: errorFlag,iRegLow,iRegUpp
  CHARACTER(CHRLEN) :: RCSIdentString
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CreateRegions.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction(global,'RFLU_CreateRegions',&
  'RFLU_CreateRegions.F90')

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  iRegUpp = global%nRegionsLocal
  
  IF ( iLev == 1 ) THEN ! finest level        
    iRegLow = 0
  ELSE ! coarse level
    iRegLow = 1
  END IF ! iLev

  ALLOCATE(levels(iLev)%regions(iRegLow:iRegUpp),STAT=errorFlag)
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'levels%regions')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************
      
  CALL DeregisterFunction(global)  
    
END SUBROUTINE RFLU_CreateRegions

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CreateRegions.F90,v $
! Revision 1.3  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/06/09 14:04:05  haselbac
! Initial revision
!
!******************************************************************************








