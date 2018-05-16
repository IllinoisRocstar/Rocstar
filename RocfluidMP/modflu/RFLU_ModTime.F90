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
! Purpose: Suite of routines to manage time.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModTime.F90,v 1.4 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTime

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_GetTimeRK, & 
            RFLU_SetTimeRK
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModTime.F90,v $ $Revision: 1.4 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
 
 
 



! ******************************************************************************
!
! Purpose: Get time within RK stages.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  FUNCTION RFLU_GetTimeRK(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Function   
! ==============================================================================

    REAL(RFREAL) :: RFLU_GetTimeRK

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    RFLU_GetTimeRK = global%currentTimeRK
     
! ******************************************************************************
!   End
! ******************************************************************************

  END FUNCTION RFLU_GetTimeRK





! ******************************************************************************
!
! Purpose: Set time within RK stages.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iStage      Index of RK stage
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_SetTimeRK(pRegion,iStage)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    INTEGER, INTENT(IN) :: iStage
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetTimeRK',&
  'RFLU_ModTime.F90')  

! ******************************************************************************
!   
! ******************************************************************************

    global%currentTimeRK = & 
      global%currentTime + global%dtMin*(pRegion%mixtInput%trk(iStage) & 
                                       - pRegion%mixtInput%trk(     1))

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SetTimeRK

 
 
 
 
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModTime


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTime.F90,v $
! Revision 1.4  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2005/03/31 17:22:54  haselbac
! Initial revision
!
! ******************************************************************************







