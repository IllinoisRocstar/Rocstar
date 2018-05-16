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
! Purpose: Set residuals to zero.
!
! Description: None.
!
! Input:
!   region        Data of current region.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: ZeroResidualsMP.F90,v 1.3 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ZeroResidualsMP(region)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), TARGET :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'ZeroResidualsMP',&
  'ZeroResidualsMP.F90')

! ******************************************************************************
! Zero residuals
! ******************************************************************************

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    DO iVar = 1,region%specInput%nSpecies
      IF ( region%specInput%specType(iVar)%frozenFlag .EQV. .TRUE. ) THEN
        region%spec%rhs(iVar,:) = 0.0_RFREAL
      END IF ! region%specInput%specType
    END DO ! iVar
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ZeroResidualsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ZeroResidualsMP.F90,v $
! Revision 1.3  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:35  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/12/01 00:06:31  wasistho
! added StatBuildVersionString
!
! Revision 1.2  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************







