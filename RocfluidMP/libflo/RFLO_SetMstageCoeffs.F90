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
! Purpose: set coefficients of an explicit multistage scheme.
!
! Description: none.
!
! Input: input = type of space and time discretization.
!
! Output: input    = stage and blending coefficients, dissipation evaluation
!         nrkSteps = number of stages.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_SetMstageCoeffs.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SetMstageCoeffs( global,input,nrkSteps )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY : t_mixt_input
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: nrkSteps

  TYPE(t_mixt_input) :: input

  TYPE(t_global), POINTER :: global

!******************************************************************************

  IF (input%timeScheme == TST_HYB5RK) THEN
    nrkSteps = 5
    IF (input%spaceDiscr == DISCR_CEN_SCAL) THEN   ! central scheme
      input%ark(1)   = 0.2500_RFREAL
      input%ark(2)   = 0.1667_RFREAL
      input%ark(3)   = 0.3750_RFREAL
      input%ark(4)   = 0.5000_RFREAL
      input%ark(5)   = 1.0000_RFREAL
      input%betrk(1) = 1.00_RFREAL
      input%betrk(2) = 0.00_RFREAL
      input%betrk(3) = 0.56_RFREAL
      input%betrk(4) = 0.00_RFREAL
      input%betrk(5) = 0.44_RFREAL
      input%ldiss(1) = 1
      input%ldiss(2) = 0
      input%ldiss(3) = 1
      input%ldiss(4) = 0
      input%ldiss(5) = 1
    ELSE                                           ! upwind scheme
      IF (input%flowModel == FLOW_EULER) THEN      ! inviscid flow
        IF (input%spaceOrder == DISCR_ORDER_1) THEN 
          input%ark(1)   = 0.0533_RFREAL
          input%ark(2)   = 0.1263_RFREAL
          input%ark(3)   = 0.2375_RFREAL
          input%ark(4)   = 0.4414_RFREAL
          input%ark(5)   = 1.0000_RFREAL
          input%betrk(1) = 1.00_RFREAL
          input%betrk(2) = 1.00_RFREAL
          input%betrk(3) = 1.00_RFREAL
          input%betrk(4) = 1.00_RFREAL
          input%betrk(5) = 1.00_RFREAL
          input%ldiss(1) = 1
          input%ldiss(2) = 1
          input%ldiss(3) = 1
          input%ldiss(4) = 1
          input%ldiss(5) = 1
        ELSE
          input%ark(1)   = 0.0695_RFREAL
          input%ark(2)   = 0.1602_RFREAL
          input%ark(3)   = 0.2898_RFREAL
          input%ark(4)   = 0.5060_RFREAL
          input%ark(5)   = 1.0000_RFREAL
          input%betrk(1) = 1.00_RFREAL
          input%betrk(2) = 1.00_RFREAL
          input%betrk(3) = 1.00_RFREAL
          input%betrk(4) = 1.00_RFREAL
          input%betrk(5) = 1.00_RFREAL
          input%ldiss(1) = 1
          input%ldiss(2) = 1
          input%ldiss(3) = 1
          input%ldiss(4) = 1
          input%ldiss(5) = 1
        ENDIF
      ELSE                                         ! viscous flow
        input%ark(1)   = 0.2742_RFREAL
        input%ark(2)   = 0.2067_RFREAL
        input%ark(3)   = 0.5020_RFREAL
        input%ark(4)   = 0.5142_RFREAL
        input%ark(5)   = 1.0000_RFREAL
        input%betrk(1) = 1.00_RFREAL
        input%betrk(2) = 0.00_RFREAL
        input%betrk(3) = 0.56_RFREAL
        input%betrk(4) = 0.00_RFREAL
        input%betrk(5) = 0.44_RFREAL
        input%ldiss(1) = 1
        input%ldiss(2) = 0
        input%ldiss(3) = 1
        input%ldiss(4) = 0
        input%ldiss(5) = 1
      ENDIF
    ENDIF

  ELSE                                ! unsteady, classical Runge-Kutta
    SELECT CASE ( global%rkScheme ) 
      CASE ( RK_SCHEME_4_CLASSICAL ) 
        nrkSteps       = 4
        input%ark(1)   = 0.5_RFREAL
        input%ark(2)   = 0.5_RFREAL
        input%ark(3)   = 1.0_RFREAL
        input%ark(4)   = 1._RFREAL/6._RFREAL
        input%betrk(:) = 1.0_RFREAL
        input%grk(1)   = 1.0_RFREAL
        input%grk(2)   = 2.0_RFREAL
        input%grk(3)   = 2.0_RFREAL
        input%grk(4)   = 1.0_RFREAL
        input%trk(1)   = 0.5_RFREAL
        input%trk(2)   = 0.5_RFREAL
        input%trk(3)   = 1.0_RFREAL
        input%trk(4)   = 1.0_RFREAL
        input%ldiss(:) = 1
        input%cfl      = MIN(input%cfl,3._RFREAL)   ! stability margin
        input%smoocf   = -1._RFREAL                 ! no residual smoothing
      CASE ( RK_SCHEME_3_WRAY ) 
        nrkSteps       = 3
        input%ark(1)   = 8.0_RFREAL/15.0_RFREAL
        input%ark(2)   = 5.0_RFREAL/12.0_RFREAL
        input%ark(3)   = 0.75_RFREAL
        input%betrk(:) = 0.0_RFREAL
        input%grk(1)   =  0.0_RFREAL
        input%grk(2)   = 17.0_RFREAL/25.0_RFREAL
        input%grk(3)   =  5.0_RFREAL/ 9.0_RFREAL        
        input%trk(1)   = 8.0_RFREAL/15.0_RFREAL 
        input%trk(2)   = 2.0_RFREAL/3.0_RFREAL 
        input%trk(3)   = 1.0_RFREAL 
        input%ldiss(:) = 1
        input%cfl = MIN(input%cfl,2.0_RFREAL)       ! stability margin
        input%smoocf   = -1._RFREAL                 ! no residual smoothing                 
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%rkScheme
  ENDIF

END SUBROUTINE RFLO_SetMstageCoeffs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SetMstageCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.2  2004/11/17 16:22:43  haselbac
! Added global as input and setting of RK3 coeffs
!
! Revision 1.1  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
!******************************************************************************






