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
! Purpose: read in user input related to random number generation
!
! Description: none.
!
! Input: user input file.
!
! Output: global%randSeedOffset
!
! Notes:
!
! * The value of global%randSeedOffset should be 0 typically (the default).
!   It can be set to 1, 2, 3, etc. to perform runs in the same configuration,
!   but with the random number generators (in each region) seeded differently.
!
! * Example input section:
!
!-----
! # RANDOM
! SEED_OFFSET  0  ! Offset for seed of RNG (default = 0, otherwise: 1,2,3,etc)
! #
!-----
!   
!******************************************************************************
!
! $Id: ReadRandomSection.F90,v 1.5 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadRandomSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX = 10

  CHARACTER(CHRLEN)     :: keys(NKEYS_MAX)

  INTEGER :: iKeySeedOffset,iKeySeedType,nKeys

  LOGICAL :: defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadRandomSection',&
  'ReadRandomSection.F90' )

! begin -----------------------------------------------------------------------

! define keys

  iKeySeedOffset = 1
  iKeySeedType   = 2
  nKeys = 2

  keys(iKeySeedOffset) = 'SEEDOFFSET'
  keys(iKeySeedType)   = 'SEEDTYPE'
  
  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  CALL ReadSection( global,IF_INPUT,nKeys,keys,vals,defined )

  IF (defined(iKeySeedOffset).eqv..true.) global%randSeedOffset = &
    ABS(NINT(vals(iKeySeedOffset)))

  IF (defined(iKeySeedType).eqv..true.) global%randSeedType = &
    ABS(NINT(vals(iKeySeedType)))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ReadRandomSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadRandomSection.F90,v $
! Revision 1.5  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.2  2005/12/01 17:08:07  fnajjar
! Added reading of randSeedType and code cleanup to remove underscore from OFFSET
!
! Revision 1.1  2004/12/01 16:50:46  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
!******************************************************************************







