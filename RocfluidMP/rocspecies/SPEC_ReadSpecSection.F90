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
! ******************************************************************************
!
! Purpose: Read user input applicable to all species.
!
! Description: None.
!
! Input:
!   regions                Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_ReadSpecSection.F90,v 1.5 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_ReadSpecSection(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

#ifdef RFLU
  USE ModInterfaces, ONLY: ReadSection
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NKEYS = 2

  CHARACTER(CHRLEN) :: keys(NKEYS)
  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: defined(NKEYS),usedFlag
  INTEGER :: errorFlag
  REAL(RFREAL) :: vals(NKEYS)
  TYPE(t_global), POINTER :: global

#ifdef RFLU
  INTEGER :: iReg
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_ReadSpecSection.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_ReadSpecSection',&
  'SPEC_ReadSpecSection.F90')

! ******************************************************************************
! Read user input for species
! ******************************************************************************

  keys(1) = 'USED'
  keys(2) = 'NSPECIES'

#ifdef RFLU
  CALL ReadSection(regions(1)%global,IF_INPUT,NKEYS,keys,vals,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN
    usedFlag = (NINT(vals(1)) /= 0)   ! if keyword USED appears, evaluate
  ELSE
    usedFlag = .TRUE.                 ! if it does not appear, default = true
  END IF ! defined(1)

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    regions(iReg)%specInput%usedFlag = usedFlag
  END DO ! iReg

  IF ( (defined(2) .EQV. .TRUE.) .AND. (usedFlag .EQV. .TRUE.) ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%specInput%nSpecies = NINT(vals(2))
      regions(iReg)%global%nSpecies    = regions(iReg)%specInput%nSpecies
    END DO ! iReg
  END IF ! defined
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_ReadSpecSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ReadSpecSection.F90,v $
! Revision 1.5  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/10/05 20:10:38  haselbac
! Added setting of global%nSpecies
!
! Revision 1.2  2004/07/28 15:31:34  jferry
! added USED field to SPECIES input section
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







