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
! Purpose: update values for Whitenoise TBC
!
! Description: none.
!
! Input: pointer to TBC
!
! Output: modifies TBC data
!
! Notes:
!
! * Example input section:
!
!-----
! # TBC_WHITENOISE
! NOSLIP    TWALL  ! BC and variable to which TBC applies
! BLOCK     0  0   ! applies to block ... (0 0 = to all)
! PATCH     0  0   ! applies to patch ... (0 0 = to all patches from blocks)
! ONTIME    5.e-5  ! time to start using this TBC
! OFFTIME   1.e-2  ! time to stop  using this TBC
! AMP       0.2    ! maximum amplitude of added stochastic component
! SUBSTEP   1      ! 0 = new value each step (default), 1 = each sub-step
! #
!-----
!
! * The value used in the boundary condition is the one input times the factor
!   1 + AMP*xi (where xi is a uniformly generated random number on [-1,1]),
!   provided OFFTIME <= t <= ONTIME.  The factor is different for every cell.
!
! * The values produced are discontinuous in time (for value continuous in time,
!   use the STOCHASTIC TBC instead).  They are independent in space.
!
!******************************************************************************
!
! $Id: UpdateTbcWhitenoise.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateTbcWhitenoise( region,tbc )

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

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'UpdateTbcWhitenoise',&
  'UpdateTbcWhitenoise.F90' )

  CALL RandUniform(tbc%bvals(TBCSTO_VAL,:),region%randData)

  tbc%bvals(TBCSTO_VAL,:) = tbc%params(TBCDAT_AMP) * &
    (2.0_RFREAL*tbc%bvals(TBCSTO_VAL,:) -  1.0_RFREAL)

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE UpdateTbcWhitenoise

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateTbcWhitenoise.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:52:14  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
! Revision 1.8  2003/06/10 22:54:23  jferry
! Added documentation for input section
!
! Revision 1.7  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.6  2003/02/17 20:53:11  jferry
! corrected misplaced assignment of global
!
! Revision 1.5  2003/02/17 19:31:12  jferry
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







