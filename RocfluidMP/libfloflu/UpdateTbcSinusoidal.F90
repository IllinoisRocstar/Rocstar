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
! Purpose: update values for Sinusoidal TBC
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
! # TBC_SINUSOIDAL
! OUTFLOW   PRESS  ! BC and variable to which TBC applies
! BLOCK     0  0   ! applies to block ... (0 0 = to all)
! PATCH     0  0   ! applies to patch ... (0 0 = to all patches from blocks)
! ONTIME    1.e-3  ! time to start using this TBC
! OFFTIME   2.e-3  ! time to stop  using this TBC
! AMP       0.2    ! amplitude of sinusoid
! FREQ      1.e4   ! frequency of sinusoid
! PHASE     0.0    ! argument of sin() for t=0
! #
!-----
!
! * The value used in the boundary condition is the one input times the factor
!   (1 + AMP*sin(FREQ*t+PHASE)), provided OFFTIME <= t <= ONTIME.
!
!******************************************************************************
!
! $Id: UpdateTbcSinusoidal.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateTbcSinusoidal( global,tbc,t )

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

!******************************************************************************

  CALL RegisterFunction( global,'UpdateTbcSinusoidal',&
  'UpdateTbcSinusoidal.F90' )

  tbc%svals(TBCSTO_VAL) = tbc%params(TBCDAT_AMP) * &
    sin(tbc%params(TBCDAT_FREQ) * t + tbc%params(TBCDAT_PHASE))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE UpdateTbcSinusoidal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateTbcSinusoidal.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:10  haselbac
! Initial revision after changing case
!
! Revision 1.6  2003/06/10 22:52:48  jferry
! Added documentation for input section
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
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
! Revision 1.1  2002/09/17 13:42:59  jferry
! Added Time-dependent boundary conditions
!
!******************************************************************************







