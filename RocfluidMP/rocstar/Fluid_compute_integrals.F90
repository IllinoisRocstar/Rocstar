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
! Purpose: Compute integrals for GENx checking.
!
! Description: None.
!
! Input:
!    globalGenx		global data structure
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: Fluid_compute_integrals.F90,v 1.4 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Fluid_compute_integrals(globalGenx,integ)

  USE ModRocstar, ONLY: t_globalGenx

#ifdef RFLO
  USE ModInterfaces, ONLY: RFLO_ComputeIntegralValues
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: RFLU_ComputeIntegralValues
#endif
  
  IMPLICIT NONE

  INCLUDE 'rocmanf90.h'

! ... parameters
  DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
  TYPE(t_globalGenx), POINTER :: globalGenx

!******************************************************************************

#ifdef RFLO
  CALL RFLO_ComputeIntegralValues(globalGenx%levels(1)%regions,integ)
#endif
#ifdef RFLU
  CALL RFLU_ComputeIntegralValues(globalGenx%levels(1)%regions,integ)
#endif


END SUBROUTINE Fluid_compute_integrals

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Fluid_compute_integrals.F90,v $
! Revision 1.4  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/02/26 04:06:02  wasistho
! added RFLO_ComputeIntegralValues
!
! Revision 1.1  2004/12/01 21:23:40  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/11/15 21:27:13  haselbac
! Initial revision
!
!******************************************************************************






