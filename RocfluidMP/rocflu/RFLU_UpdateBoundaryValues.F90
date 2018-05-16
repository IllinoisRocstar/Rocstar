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
! Purpose: Update boundary values.
!
! Description: None.
!
! Input: 
!   region     Data of current region
!   istage     Current RK stage
!
! Output: 
!   None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_UpdateBoundaryValues.F90,v 1.6 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_UpdateBoundaryValues(region,istage)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetBoundaryValues, &
                           RFLU_PutBoundaryValuesAlpha, & 
                           UpdateTbc
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: istage
  TYPE(t_region) :: region  

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: subdt,time
  REAL(RFREAL) :: trk(5)
  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: alpha
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_UpdateBoundaryValues.F90,v $ $Revision: 1.6 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_UpdateBoundaryValues',&
  'RFLU_UpdateBoundaryValues.F90')

  trk(:) = region%mixtInput%trk(:)

! ******************************************************************************
! Compute time
! ******************************************************************************

  time = global%currentTime + global%dtMin*trk(istage)

#ifdef GENX 
! ******************************************************************************
! Fill incoming buffers
! ******************************************************************************
   
  alpha = (time-global%timeStamp)/global%dTimeSystem

  CALL COM_call_function(global%genxHandleBc,2,alpha,1)
  CALL RFLU_PutBoundaryValuesAlpha(region)
  CALL COM_call_function(global%genxHandleBc,2,alpha,2)
  CALL RFLU_GetBoundaryValues(region)
#endif

! ******************************************************************************
! Fill in time-dependent boundary condition data 
! ******************************************************************************

  IF ( istage > 1 ) THEN 
    subdt = global%dtMin*(trk(istage) - trk(istage-1))
  ELSE
    subdt = 0.0_RFREAL
  END IF ! istage

  CALL UpdateTbc(region,time,subdt,istage==global%nrkSteps)  

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_UpdateBoundaryValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_UpdateBoundaryValues.F90,v $
! Revision 1.6  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/04/15 15:07:28  haselbac
! Removed Charm/FEM stuff, cosmetics
!
! Revision 1.3  2004/01/22 16:04:34  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.2  2003/11/25 21:07:44  haselbac
! Temporarily disable UpdateTbc - leads to core dump
!
! Revision 1.1  2003/10/03 20:48:52  haselbac
! Initial revision
!
! ******************************************************************************







