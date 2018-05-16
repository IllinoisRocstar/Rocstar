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
! Purpose: Read user input, store it in the data structure and check.
!
! Description: None.
!
! Input:
!   regions                Data associated with regions
!   iSpecType                Species type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_ReadSpecTypeSection.F90,v 1.10 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_ReadSpecTypeSection(regions,iSpecType)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModSpecies, ONLY: t_spec_type

#ifdef RFLU
  USE ModInterfaces, ONLY: ReadBothSection, &
                           ReadSection
  USE ModInterfacesInteract, ONLY: INRT_SetMaterial
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iSpecType
  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NKEYS = 8, NSTRKEYS = 1

  CHARACTER(CHRLEN) :: keys(NKEYS),strKeys(NSTRKEYS),strVals(NSTRKEYS)
  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: defined(NKEYS),strDefined(NSTRKEYS)
  REAL(RFREAL) :: vals(NKEYS)
  TYPE(t_global), POINTER :: global
  TYPE(t_spec_type), POINTER :: pSpecType

#ifdef RFLU
  INTEGER :: iReg
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_ReadSpecTypeSection.F90,v $ $Revision: 1.10 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_ReadSpecTypeSection',&
  'SPEC_ReadSpecTypeSection.F90')

! ******************************************************************************
! Read user input for species
! ******************************************************************************

  keys(1) = 'FROZENFLAG'
  keys(2) = 'INITVAL'
  keys(3) = 'SOURCETYPE'
  keys(4) = 'SCHMIDTNO'
  keys(5) = 'DIAMETER'
  keys(6) = 'PUFFFACTOR'
  keys(7) = 'VELOCITYMETHOD'
  keys(8) = 'SETTLINGFLAG'

  strKeys(1) = 'MATERIAL'

#ifdef RFLU
! ==============================================================================
! Read SPEC_TYPE section
! ==============================================================================

  CALL ReadBothSection(global,IF_INPUT,NKEYS,NSTRKEYS,keys,strKeys,vals, &
                       strVals,defined,strDefined)


! ==============================================================================
! Set variables
! ==============================================================================

! ------------------------------------------------------------------------------
! Frozen flag
! ------------------------------------------------------------------------------

  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        IF ( NINT(vals(1)) == 0 ) THEN
          regions(iReg)%specInput%specType(iSpecType)%frozenFlag = .FALSE.
        ELSE
          regions(iReg)%specInput%specType(iSpecType)%frozenFlag = .TRUE.
        END IF ! NINT
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Initial value
! ------------------------------------------------------------------------------

  IF ( defined(2) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%initVal = vals(2)
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Source type
! ------------------------------------------------------------------------------

  IF ( defined(3) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%sourceType = NINT(vals(3))
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Schmidt number
! ------------------------------------------------------------------------------

  IF ( defined(4) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%schmidtNumber = vals(4)
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Diameter
! ------------------------------------------------------------------------------

  IF ( defined(5) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%diameter = vals(5)
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Puff factor
! ------------------------------------------------------------------------------

  IF ( defined(6) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%puffFactor = vals(6)
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined

! ------------------------------------------------------------------------------
! Velocity method
! ------------------------------------------------------------------------------

  IF ( defined(7) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        regions(iReg)%specInput%specType(iSpecType)%velocityMethod = &
          NINT(vals(7))
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined
  
! ------------------------------------------------------------------------------
! Settling flag
! ------------------------------------------------------------------------------

  IF ( defined(8) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        IF ( NINT(vals(8)) == 0 ) THEN
          regions(iReg)%specInput%specType(iSpecType)%settlingFlag = .FALSE.
        ELSE
          regions(iReg)%specInput%specType(iSpecType)%settlingFlag = .TRUE.
        END IF ! NINT
      END IF ! iSpecType
    END DO ! iReg
  END IF ! defined 

! ------------------------------------------------------------------------------
! Set material
! ------------------------------------------------------------------------------

  IF ( strDefined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      IF ( iSpecType <= regions(iReg)%specInput%nSpecies ) THEN
        pSpecType => regions(iReg)%specInput%specType(iSpecType)

        CALL INRT_SetMaterial(global,pSpecType%pMaterial,strVals(1))
      END IF ! iSpecType
    END DO ! iReg
  END IF ! strDefined
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_ReadSpecTypeSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ReadSpecTypeSection.F90,v $
! Revision 1.10  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2005/11/10 02:39:01  haselbac
! Added support for settling flag
!
! Revision 1.7  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.6  2004/07/28 15:31:34  jferry
! added USED field to SPECIES input section
!
! Revision 1.5  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/04/19 20:21:54  haselbac
! Added code to read source type
!
! Revision 1.3  2004/02/02 22:42:27  haselbac
! Added reading of material name and setting of material pointer
!
! Revision 1.2  2004/01/29 23:00:44  haselbac
! Added reading of Schmidt number
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







