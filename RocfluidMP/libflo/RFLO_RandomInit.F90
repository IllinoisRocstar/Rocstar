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
! Notes:
!
!   Uses positive seeds.
!
!   The seed 0 is forbidden; negative seeds are used for unstructured regions.
!
!******************************************************************************
!
! $Id: RFLO_RandomInit.F90,v 1.4 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_RandomInit( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModRandom,     ONLY : RandomSeed
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: iclock1,seed,seedClock
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_RandomInit',&
  'RFLO_RandomInit.F90' )

! ******************************************************************************
! Set seed type
! ******************************************************************************

  SELECT CASE(global%randSeedType)
    CASE(RAND_SEED_TYPE_FIXED)
      seedClock = 0
    CASE(RAND_SEED_TYPE_CLOCK)
      CALL system_clock(iclock1)
      seedClock = iclock1
  END SELECT ! global%randSeedType

! ******************************************************************************
! Loop over regions
! ******************************************************************************

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      seed = iReg + global%nRegions * global%randSeedOffset +seedClock
      IF (seed <= 0) THEN  ! RFLO seeds must be positive
        CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__ )
      ENDIF
      CALL RandomSeed(seed,regions(iReg)%randData)
    ENDIF ! procid
  ENDDO

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_RandomInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_RandomInit.F90,v $
! Revision 1.4  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/12/01 17:25:04  fnajjar
! Added capability for clock-based random seed
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.3  2003/11/21 22:30:04  fnajjar
! Added global seed offset for Random Number Generator
!
! Revision 1.2  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.1  2003/02/17 19:31:11  jferry
! Implemented portable random number generator ModRandom
!
!******************************************************************************







