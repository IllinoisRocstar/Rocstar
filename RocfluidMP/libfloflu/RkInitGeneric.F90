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
! Purpose: Initialize Runge-Kutta scheme.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   cv          Conserved variables
!   cvOld       Old conserved variables
!   diss        Dissipation residual
!
! Output: 
!   cvOld       Old conserved variables
!   diss        Dissipation residual
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkInitGeneric.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkInitGeneric(region,iStage,icBeg,icEnd,ivBeg,ivEnd,cv,cvOld,diss)

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

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,diss
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

  RCSIdentString = '$RCSfile: RkInitGeneric.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction(global,'RkInitGeneric',&
  'RkInitGeneric.F90')

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            cvOld(iv,ic) = cv(iv,ic)
            diss(iv,ic)  = 0.0_RFREAL
          END DO ! iv
        END DO ! ic
      ELSE
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            diss(iv,ic) = 0.0_RFREAL
          END DO ! iv      
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          cvOld(iv,ic) = cv(iv,ic)
          diss(iv,ic)  = 0.0_RFREAL
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RkInitGeneric

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitGeneric.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:51:02  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/11/17 16:23:40  haselbac
! Added initialization for RK3
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************







