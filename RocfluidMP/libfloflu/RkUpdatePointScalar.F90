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
! Purpose: Updates particle velocity with classical 4-stage Runge-Kutta method.
!
! Description: None.
!
! Input: 
!   region	Region data
!   iStage	Runge-Kutta stage
!   ivBeg	Beginning index for variable update
!   ivEnd	Ending index for variable update
!   var		Conserved variables
!   varOld	Old conserved variables
!   rhs		Residual
!   rhsSum	Residual sum
!
!
! Output: 
!   var		Variables
!   rhsSum	Residual sum
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkUpdatePointScalar.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkUpdatePointScalar(region,iStage,ivBeg,ivEnd,var,varOld,rhs,rhsSum)

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

  TYPE(t_region) :: region
  INTEGER, INTENT(IN) :: iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:), POINTER :: var,varOld,rhs,rhsSum

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iv
  REAL(RFREAL) :: fac
  REAL(RFREAL) :: ark(5),grk(5)
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkUpdatePointScalar.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction(global,'RkUpdatePointScalar',&
  'RkUpdatePointScalar.F90')

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! *****************************************************************************
! Update
! *****************************************************************************

  fac = ark(iStage)*global%dtMin

  SELECT CASE ( global%rkScheme ) 
  CASE ( RK_SCHEME_4_CLASSICAL ) 
    IF ( iStage == 1 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE IF ( iStage == global%nrkSteps ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) + rhsSum(iv))
      END DO ! iv
    ELSE
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhsSum(iv) + grk(iStage)*rhs(iv)
      END DO ! iv
    END IF ! iStage
  CASE ( RK_SCHEME_3_WRAY ) 
    IF ( iStage == 1 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE IF ( iStage == 2 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) - grk(iStage)*rhsSum(iv))
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) - grk(iStage)*rhsSum(iv))
      END DO ! iv
    END IF ! iStage      
  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

 CALL DeregisterFunction(global)

END SUBROUTINE RkUpdatePointScalar

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdatePointScalar.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/19 15:37:27  mparmar
! Initial revision
!
! ******************************************************************************







