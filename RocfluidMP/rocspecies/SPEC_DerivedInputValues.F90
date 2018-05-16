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
! Purpose: Set derived input values for species.
!
! Description: None.
!
! Input:
!   regions     Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_DerivedInputValues.F90,v 1.12 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_DerivedInputValues(regions)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region  
  USE ModSpecies, ONLY: t_spec_type

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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg,iSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_spec_type), POINTER :: pSpecType

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_DerivedInputValues.F90,v $ $Revision: 1.12 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_DerivedInputValues',&
  'SPEC_DerivedInputValues.F90')

! ******************************************************************************
! Loop over regions
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    regions(iReg)%spec%nSpecEqs = regions(iReg)%specInput%nSpecies

! ==============================================================================
!   Loop over species
! ==============================================================================

    regions(iReg)%specInput%nSpeciesEE = 0

    DO iSpec = 1,regions(iReg)%specInput%nSpecies
      pSpecType => regions(iReg)%specInput%specType(iSpec)

! ------------------------------------------------------------------------------
!     Set source flag
! ------------------------------------------------------------------------------

      SELECT CASE ( pSpecType%sourceType )
        CASE ( SPEC_SOURCE_TYPE_NONE )
          regions(iReg)%specInput%sourceFlag = .FALSE.  
        CASE ( SPEC_SOURCE_TYPE_CHEM )
          regions(iReg)%specInput%sourceFlag = .TRUE.
        CASE ( SPEC_SOURCE_TYPE_CAVI )
          regions(iReg)%specInput%sourceFlag = .TRUE.
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__, &
                         'Invalid input value for SOURCETYPE')
      END SELECT ! sourceType

! ------------------------------------------------------------------------------
!     Set discreteFlag and related variables
! ------------------------------------------------------------------------------

      pSpecType%discreteFlag = (pSpecType%pMaterial%phase /= 1) ! i.e., not Gas

      IF ( pSpecType%discreteFlag .EQV. .FALSE. ) THEN 
        pSpecType%iCont = 1
      ELSE 
        pSpecType%iCont = 0
      END IF ! pSpecType%discreteFlag

! --- Discrete species ---------------------------------------------------------

      IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
        pSpecType%effectiveDensity = pSpecType%pMaterial%dens/  &
                                     pSpecType%puffFactor
        pSpecType%effectiveVolume  = global%pi*pSpecType%diameter**3/ &
                                     6.0_RFREAL
        pSpecType%materialVolume   = pSpecType%effectiveVolume/ &
                                     pSpecType%puffFactor

! ----- Include pressure correction: (1 - beta) factor. We have chosen *not* 
!       to set this to 0 for the fluid vel case

        pSpecType%tauCoefficient = pSpecType%diameter**2* &
          (pSpecType%effectiveDensity-1.0_RFREAL)/18.0_RFREAL

! ----- Do not waste Eq Eul method on cases with tau *exactly* equal to 0

        IF ( pSpecType%tauCoefficient == 0.0_RFREAL ) THEN
          pSpecType%velocityMethod = SPEC_METHV_FLUIDVEL
        END IF ! pSpecType%tauCoefficient

! --- Continuous species -------------------------------------------------------

      ELSE
        pSpecType%diameter         = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pSpecType%puffFactor       = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pSpecType%effectiveDensity = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pSpecType%effectiveVolume  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pSpecType%materialVolume   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pSpecType%tauCoefficient   = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

        IF ( pSpecType%velocityMethod /= SPEC_METHV_FLUIDVEL ) THEN
          CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
            'Continuum species must use fluid velocity for advection')
        END IF ! velocityMethod
      END IF ! discreteFlag

! ------------------------------------------------------------------------------
!     Set indSd to address substantial derivative
! ------------------------------------------------------------------------------

      SELECT CASE ( pSpecType%velocityMethod )
        CASE ( SPEC_METHV_FLUIDVEL )
          regions(iReg)%mixtInput%indSd = 0     
        CASE ( SPEC_METHV_EQEUL )
          regions(iReg)%mixtInput%indSd = 1

          regions(iReg)%specInput%nSpeciesEE = & 
            regions(iReg)%specInput%nSpeciesEE + 1
          pSpecType%iSpec2iSpecEEv = regions(iReg)%specInput%nSpeciesEE
          pSpecType%iSpecEEv2iSpec = iSpec         
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__, &
                         'Invalid input value for VELOCITYMETHOD')
      END SELECT ! velocityMethod
    END DO ! iSpec

! - require viscosity to be computed if Equilibrium Eulerian velocity is used

    IF ( regions(iReg)%mixtInput%indSd == 1 ) THEN
      regions(iReg)%mixtInput%computeTv = .TRUE.
      
      IF ( regions(iReg)%mixtInput%nTv < 2 ) THEN
        regions(iReg)%mixtInput%nTv = 2
      END IF ! nTv
    END IF ! indSd
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_DerivedInputValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_DerivedInputValues.F90,v $
! Revision 1.12  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.9  2006/03/30 20:50:41  haselbac
! Added setting of sourceFlag for cavitation src term
!
! Revision 1.8  2005/11/27 01:53:40  haselbac
! Added setting of EEv-related variables
!
! Revision 1.7  2005/11/14 17:02:17  haselbac
! Clean-up, added setting of iCont
!
! Revision 1.6  2004/08/02 16:38:40  jferry
! Added pressure correction to tau coefficient
!
! Revision 1.5  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.4  2004/07/28 15:31:34  jferry
! added USED field to SPECIES input section
!
! Revision 1.3  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.2  2004/04/01 21:31:18  haselbac
! Added setting of sourceFlag
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







