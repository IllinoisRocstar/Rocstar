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
! Purpose: Initialize array sd, which holds the substantial derivative used
!          in computing Equilibrium Eulerian velocities
!
! Description: None.
!
! Input:
!   region       Region data
!   icBeg        Beginning index for cell update
!   icEnd        Ending index for cell update
!   ivBeg        Beginning index for variable update
!   ivEnd        Ending index for variable update
!
! Output:
!   sd           Substantial derivative
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkInitSD.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkInitSD(region,icBeg,icEnd,ivBeg,ivEnd,sd)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icBeg,icEnd,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:,:), POINTER :: sd
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkInitSD.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction(global,'RkInitSD',&
  'RkInitSD.F90')

! *****************************************************************************
! Initialize substantial derivative
! *****************************************************************************

  DO ic = icBeg,icEnd
    DO iv = ivBeg,ivEnd
      sd(iv,ic) = 0.0_RFREAL
    END DO ! iv
  END DO ! ic

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RkInitSD

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitSD.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:51:06  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
!******************************************************************************







