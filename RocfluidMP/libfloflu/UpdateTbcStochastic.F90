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
! Purpose: update values for Stochastic TBC
!
! Description: none.
!
! Input: pointer to TBC, substep dt
!
! Output: modifies TBC data
!
! Notes:
!
! * Example input section:
!
!-----
! # TBC_STOCHASTIC
! PEUL_INJECT  FRAC2 ! BC and variable to which TBC applies
! BLOCK     0  0   ! applies to block ... (0 0 = to all)
! PATCH     0  0   ! applies to patch ... (0 0 = to all patches from blocks)
! ONTIME    5.e-5  ! time to start using this TBC
! OFFTIME   1.e-2  ! time to stop  using this TBC
! AMP       0.2    ! standard deviation of log of the stochastic factor
! TIMECOR   1.e-4  ! integral time scale of autocorrelation
! SHAPE     1.0    ! shape:  recommended value = 1.0 is default
! MINCUT    0.2    ! minimal value for stochastic factor (< 0 for none)
! MAXCUT    5.0    ! maximal value for stochastic factor (< 0 for none)
! #
!-----
!
! * The value used in the boundary condition is the one input times the factor
!   produced by this routine, provided OFFTIME <= t <= ONTIME.  The factor is
!   different for every cell, and is always positive.
!
! * The values produced are continuous and, in fact, smooth in time.  They are
!   independent in space.
!
! * TIMECOR is the decorrelation time scale of the process.
!
! * SHAPE determines what the autocorrelation function looks like.  For
!
!     0 <= SHAPE < 1, R(t) is overdamped,
!     SHAPE = 1,      R(t) is critically damped,
!     SHAPE > 1,      R(t) is underdamped (i.e., decaying but oscillatory).
!
! * For numerical reasons, the restriction 0.1 <= SHAPE <= 10 is imposed.
!
!******************************************************************************
!
! $Id: UpdateTbcStochastic.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateTbcStochastic( region,tbc,dt )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_tbcvalues
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModParameters
  USE ModRandom,     ONLY : RandUniform
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_region),    INTENT(INOUT) :: region
  TYPE(t_tbcvalues), INTENT(INOUT) :: tbc

  REAL(RFREAL),      INTENT(IN)    :: dt

! ... local variables
  REAL(RFREAL)            :: dtn,fac,e0,e1,s0,s1,ts,a0,b0,a1,b1,del
  REAL(RFREAL),   POINTER :: params(:), bvals(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'UpdateTbcStochastic',&
  'UpdateTbcStochastic.F90' )

  bvals  => tbc%bvals
  params => tbc%params

  dtn = dt/params(TBCDAT_TIMECOR)
  fac = -2.0_RFREAL/params(TBCDAT_SHAPE)

  e0 = EXP(fac*dtn*dtn)
  e1 = EXP(2.0_RFREAL*fac*dtn)
  s0 = SQRT(1.0_RFREAL - e0*e0)
  s1 = SQRT(1.0_RFREAL - e1*e1)
  ts = 0.5_RFREAL * params(TBCDAT_TIMECOR) * SQRT(params(TBCDAT_SHAPE))

  a0 = e0
  b0 = s0 * ts
  a1 = -s0*e1 / (ts*e0)
  b1 = e1/e0
  del = SQRT(12.0_RFREAL) * params(TBCDAT_AMP) * s1 / ts

  CALL RandUniform(bvals(TBCSTO_FACTOR,:),region%randData)

  bvals(TBCSTO_VAL, :)   = a0*bvals(TBCSTO_VAL,:) + b0*bvals(TBCSTO_DVAL,:)
  bvals(TBCSTO_DVAL,:)   = a1*bvals(TBCSTO_VAL,:) + b1*bvals(TBCSTO_DVAL,:) + &
                             del * (bvals(TBCSTO_FACTOR,:) - 0.5_RFREAL)
  bvals(TBCSTO_FACTOR,:) = EXP(bvals(TBCSTO_VAL,:) - &
                             0.5_RFREAL * params(TBCDAT_AMP)**2 )
  IF (params(TBCDAT_MINCUT) >= 0._RFREAL) &
    bvals(TBCSTO_FACTOR,:) = MAX(params(TBCDAT_MINCUT),bvals(TBCSTO_FACTOR,:))
  IF (params(TBCDAT_MAXCUT) >= 0._RFREAL) &
    bvals(TBCSTO_FACTOR,:) = MIN(params(TBCDAT_MAXCUT),bvals(TBCSTO_FACTOR,:))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE UpdateTbcStochastic

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateTbcStochastic.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:12  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
! Revision 1.7  2003/06/10 22:52:48  jferry
! Added documentation for input section
!
! Revision 1.6  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.5  2003/02/17 19:31:11  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.4  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.3  2002/09/25 18:29:57  jferry
! simplified TBC parameter lists
!
! Revision 1.2  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.1  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
!******************************************************************************







