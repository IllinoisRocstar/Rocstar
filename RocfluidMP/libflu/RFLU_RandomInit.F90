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
! Purpose: Initializes random number generator in each region
!
! Description: none.
!
! Input:  does not use any data already present
!
! Output: regions(:)%randData is initialized
!
! Notes: none 
!
!  Uses negative seeds.  The procedure to produce them is somewhat awkward
!  because we make no assumptions about the values regions(:)%iRegionGlobal
!
!  The seed 0 is forbidden; positive seeds are used for structured regions.

!******************************************************************************
!
! $Id: RFLU_RandomInit.F90,v 1.5 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_RandomInit( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModRandom,     ONLY : RandomSeed
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iclock1,ilo,ihi,iReg,seed,seedClock

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_RandomInit.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLU_RandomInit',&
  'RFLU_RandomInit.F90' )

! ******************************************************************************
! Set seed type
! ******************************************************************************

  SELECT CASE(global%randSeedType)
    CASE(RAND_SEED_TYPE_FIXED)
      seedClock = 0
    CASE(RAND_SEED_TYPE_CLOCK)
      CALL system_clock(iclock1)
      seedClock = -iclock1
  END SELECT ! global%randSeedType

! ******************************************************************************
! Loop over regions
!  Note: seed with negative numbers is used for RFLU
! ******************************************************************************

  ilo = LBOUND(regions,1)
  ihi = UBOUND(regions,1)

  DO iReg = ilo,ihi
    seed = iReg-ihi-1 - (ihi-ilo+1) * global%randSeedOffset + seedClock
    IF (seed >= 0) THEN  ! RFLU seeds must be negative
      CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__ )
    ENDIF ! seed

    CALL RandomSeed(seed,regions(iReg)%randData)
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_RandomInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_RandomInit.F90,v $
! Revision 1.5  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/12/01 17:09:05  fnajjar
! Added capability for clock-based random seed
!
! Revision 1.2  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.1  2003/02/17 19:31:11  jferry
! Implemented portable random number generator ModRandom
!
!
!******************************************************************************







