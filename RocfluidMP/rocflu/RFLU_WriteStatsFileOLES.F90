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
! Purpose: Write optimal LES statistics in ASCII ROCFLU format.
!
! Description: None.
!
! Input: 
!   region      region data
!
! Output: None.
!
! Notes: 
!   1. Works only for single regions right now
!
!******************************************************************************
!
! $Id: RFLU_WriteStatsFileOLES.F90,v 1.3 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteStatsFileOLES(global)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError    

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================
  
  CHARACTER(CHRLEN) :: RCSIdentString  
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  
! ******************************************************************************
! Start, set pointer
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteStatsFileOLES.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction(global,'RFLU_WriteStatsFileOLES',&
  'RFLU_WriteStatsFileOLES.F90')

! ******************************************************************************
! Write data
! ******************************************************************************

  WRITE(IF_STATS_OLES,'(6(2X,E15.8))') & 
    global%currentTime,global%enerOLES,global%dissOLES, & 
    global%uVarOLES,global%vVarOLES,global%wVarOLES

  CALL DeregisterFunction(global)
     
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_WriteStatsFileOLES


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_WriteStatsFileOLES.F90,v $
!   Revision 1.3  2008/12/06 08:44:30  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:43  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2002/09/09 16:28:02  haselbac
!   Initial revision
!
! ******************************************************************************







