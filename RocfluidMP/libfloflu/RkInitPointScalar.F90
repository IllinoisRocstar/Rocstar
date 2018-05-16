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
! Purpose: Initialize Runge-Kutta scheme particle velocity computations.
!
! Description: None.
!
! Input: 
!   region	Region data
!   iStage	Runge-Kutta stage
!   ivBeg	Beginning index for variable update
!   ivEnd	Ending index for variable update
!   var		Conserved variables
!   varOld	Old variables
!
! Output: 
!   varOld	Old variables
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkInitPointScalar.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RkInitPointScalar(region,iStage,ivBeg,ivEnd,var,varOld)

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

  INTEGER, INTENT(IN) :: iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:), POINTER :: var,varOld
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iv
  TYPE(t_global), POINTER :: global
  
! *****************************************************************************
! Start
! *****************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RkInitPointScalar',&
  'RkInitPointScalar.F90')

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
  CASE ( RK_SCHEME_4_CLASSICAL ) 
    IF ( iStage == 1 ) THEN
      DO iv = ivBeg,ivEnd
        varOld(iv) = var(iv) 
      END DO ! iv
    END IF ! iStage
  CASE ( RK_SCHEME_3_WRAY ) 
    DO iv = ivBeg,ivEnd
      varOld(iv) = var(iv) 
    END DO ! iv
  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RkInitPointScalar

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitPointScalar.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/19 15:37:30  mparmar
! Initial revision
!
! ******************************************************************************







