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
! Purpose: update values for Piecewise Constant or Piecewise Linear TBC
!
! Description: none.
!
! Input: pointer to TBC, substep time
!
! Output: modifies TBC data
!
! Notes:
!
! * Example input section:
!
!-----
! # TBC_PIECEWISE
! INJECT    MFRATE ! BC and variable to which TBC applies
! BLOCK     0  0   ! applies to block ... (0 0 = to all)
! PATCH     0  0   ! applies to patch ... (0 0 = to all patches from blocks)
! ONTIME   -1.e6   ! time to start using this TBC
! OFFTIME   1.e6   ! time to stop  using this TBC
! ORDER     0      ! 0 = piecewise constant (default), 1 = piecewise linear
! NJUMPS    4      ! number of points at which behavior changes
! #
! FRAC  0.0        ! fraction of input value of variable before first time
! TIME  0.001      ! first time at which behavior changes
! FRAC  0.1        ! next fraction attained (constant) or ramped to (linear)
! TIME  0.002      ! second time at which behavior changes
! FRAC  0.3
! TIME  0.003
! FRAC  0.6
! TIME  0.004      ! final time at which behavior changes
! FRAC  1.0        ! final value for constant case; *ignored* for linear case
! #
!-----
!
! * ONTIME and OFFTIME would not typically be used with this TBC.
!
! * ORDER = 0 is specified above, which yields
!
!   frac = 0   for         t < 0.001,
!   frac = 0.1 for 0.001 < t < 0.002,
!   frac = 0.3 for 0.002 < t < 0.003,
!   frac = 0.6 for 0.003 < t < 0.004,
!   frac = 1.0 for 0.004 < t.
!
! * If ORDER = 1 were specified instead, the above would yield
!
!   frac = 0            for         t < 0.001,
!   frac = 100.0*t-0.1  for 0.001 < t < 0.002,
!   frac = 200.0*t-0.3  for 0.002 < t < 0.003,
!   frac = 300.0*t-0.6  for 0.003 < t < 0.004,
!   frac = 0.6          for 0.004 < t.
!
!   Note that the final value (1.0) is ignored for this case.
!
! * The value used in the boundary condition is the one input times the factor
!   "frac" given above, provided OFFTIME <= t <= ONTIME.
!
!******************************************************************************
!
! $Id: UpdateTbcPiecewise.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateTbcPiecewise( global,tbc,t )

  USE ModDataTypes
  USE ModBndPatch, ONLY : t_tbcvalues
  USE ModGlobal,   ONLY : t_global
  USE ModParameters
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_global),    POINTER       :: global
  TYPE(t_tbcvalues), INTENT(INOUT) :: tbc

  REAL(RFREAL),      INTENT(IN)    :: t

! ... loop variables
  INTEGER :: i

! ... local variables
  INTEGER      :: order,njumps
  LOGICAL      :: usemaxtime
  REAL(RFREAL) :: val,val0,val1,t0,t1

!******************************************************************************

  CALL RegisterFunction( global,'UpdateTbcPiecewise',&
  'UpdateTbcPiecewise.F90')

  order  = tbc%switches(TBCSWI_ORDER)
  njumps = tbc%switches(TBCSWI_NJUMPS)

  IF (t >= tbc%params(TBCDAT_DAT0 + 2*njumps)) THEN

    usemaxtime = .true.

  ELSE IF (t < tbc%params(TBCDAT_DAT0 + 2)) THEN

    usemaxtime = .false.
    val = tbc%params(TBCDAT_DAT0 + 1) ! for both CONSTANT and LINEAR cases

  ELSE

    usemaxtime = .true.
    DO i=2,njumps
      IF (t < tbc%params(TBCDAT_DAT0 + 2*i)) THEN
        SELECT CASE (order)
        CASE (TBCOPT_CONSTANT)
          val = tbc%params(TBCDAT_DAT0 + 2*i-1) ! CONSTANT case
        CASE (TBCOPT_LINEAR)
          val0 = tbc%params(TBCDAT_DAT0 + 2*i-3)
          t0   = tbc%params(TBCDAT_DAT0 + 2*i-2)
          val1 = tbc%params(TBCDAT_DAT0 + 2*i-1)
          t1   = tbc%params(TBCDAT_DAT0 + 2*i)
          val  = ((t-t0)*val1 + (t1-t)*val0) / (t1-t0) ! LINEAR case
        END SELECT ! order
        usemaxtime = .false.
        EXIT
      ENDIF ! t
    ENDDO

  ENDIF ! tbc%params

  IF (usemaxtime) THEN

    SELECT CASE (order)
    CASE (TBCOPT_CONSTANT)
      val = tbc%params(TBCDAT_DAT0 + 2*njumps+1) ! CONSTANT case
    CASE (TBCOPT_LINEAR)
      val = tbc%params(TBCDAT_DAT0 + 2*njumps-1) ! LINEAR case
    END SELECT ! order

  ENDIF ! usemaxtime

  tbc%svals(TBCSTO_VAL) = val

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE UpdateTbcPiecewise

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateTbcPiecewise.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:07  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
!******************************************************************************







