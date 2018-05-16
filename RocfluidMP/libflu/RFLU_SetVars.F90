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
! Purpose: Wrapper for setting mixture variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetVars.F90,v 1.6 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE ModInterfaces, ONLY: RFLU_SetGasVars, & 
                           RFLU_SetDependentVars, & 
                           RFLU_SetTransportVars
                     
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVars.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVars',&
  'RFLU_SetVars.F90')

#ifdef ROCPROF
  CALL FPROFILER_BEGINS("RFLU::SetVars")
#endif

! ******************************************************************************
! Set variables
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_SetGasVars(pRegion,icgBeg,icgEnd)
  CALL RFLU_SetDependentVars(pRegion,icgBeg,icgEnd)
  
  IF ( pRegion%mixtInput%computeTv .EQV. .TRUE. ) THEN 
    CALL RFLU_SetTransportVars(pRegion,icgBeg,icgEnd)
  END IF ! pRegion%mixtInput%computeTv

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF
  CALL FPROFILER_ENDS("RFLU::SetVars")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVars.F90,v $
! Revision 1.6  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2005/07/07 22:44:24  haselbac
! Added profiling calls
!
! Revision 1.2  2005/04/15 15:06:23  haselbac
! Added range arguments, adapted calls to other routines accordingly
!
! Revision 1.1  2004/11/14 20:02:49  haselbac
! Initial revision
!
! ******************************************************************************







