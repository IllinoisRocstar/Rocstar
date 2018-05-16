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
! Purpose: Wrapper routine to get deformation.
!
! Description: None.
!
! Input: 
!   region     Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_GetDeformationWrapper.F90,v 1.4 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_GetDeformationWrapper(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetDeformation
#else
  USE ModInterfaces, ONLY: RFLU_USER_GetDeformation
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: dAlpha
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetDeformationWrapper.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_GetDeformationWrapper',&
  'RFLU_GetDeformationWrapper.F90')

! *****************************************************************************
! Get deformation
! *****************************************************************************
  
#ifdef GENX 
  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function(global%genxHandleGm,1,dAlpha)
  
  DO iReg = 1,global%nRegionsLocal
    CALL RFLU_GetDeformation(regions(iReg))
  END DO ! iReg  
#else
  DO iReg = 1,global%nRegionsLocal
    CALL RFLU_USER_GetDeformation(regions(iReg))
  END DO ! iReg
#endif     

! *****************************************************************************
! End
! *****************************************************************************
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetDeformationWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetDeformationWrapper.F90,v $
! Revision 1.4  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2003/05/09 17:01:04  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.1  2003/03/15 19:06:12  haselbac
! Initial revision
!
!******************************************************************************







