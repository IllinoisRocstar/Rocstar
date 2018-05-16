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
! Purpose: Fill GENX incoming buffers.
!
! Description: None.
!
! Input: 
!   region     Data of current region
!   istage     Current RK stage
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: RFLU_UpdateInbuffWrapper.F90,v 1.6 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_UpdateInbuffWrapper(region,istage)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetBoundaryValues, &
                           RFLU_PutBoundaryValuesAlpha
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
  REAL(RFREAL) :: time
  REAL(RFREAL) :: trk(5)
  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: alpha
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_UpdateInbuffWrapper.F90,v $ $Revision: 1.6 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_UpdateInbuffWrapper',&
  'RFLU_UpdateInbuffWrapper.F90')

#ifdef GENX 
! ******************************************************************************
! Fill incoming buffers
! ******************************************************************************
   
  trk(:) = region%mixtInput%trk(:)

  time  = global%currentTime + global%dtMin*(trk(istage) - trk(1))
  alpha = (time-global%timeStamp)/global%dTimeSystem

  CALL COM_call_function(global%genxHandleBc,2,alpha,1)
  CALL RFLU_PutBoundaryValuesAlpha(region)
  CALL COM_call_function(global%genxHandleBc,2,alpha,2)
  CALL RFLU_GetBoundaryValues(region)
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_UpdateInbuffWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_UpdateInbuffWrapper.F90,v $
! Revision 1.6  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/04/15 15:05:59  haselbac
! Removed USE RFLU_ModFEM, cosmetics
!
! Revision 1.3  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.2  2003/04/12 21:35:51  haselbac
! Clean-up
!
! Revision 1.1  2003/03/28 19:37:18  fnajjar
! Initial Import for RocfluidMP
!
! ******************************************************************************







