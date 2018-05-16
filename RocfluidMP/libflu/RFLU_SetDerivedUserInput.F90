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
! Purpose: Set values derived from user input of ALL modules.
!
! Description: None.
!
! Input:
!   regions        Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetDerivedUserInput.F90,v 1.12 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetDerivedUserInput(regions)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModMPI
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModSpecies, ONLY: t_spec_input
  USE ModDataStruct, ONLY: t_region

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_spec_input), POINTER :: pSpecInput
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

#ifdef SPEC
  INTEGER :: contCntr,iSpec
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetDerivedUserInput.F90,v $ $Revision: 1.12 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_SetDerivedUserInput',&
  'RFLU_SetDerivedUserInput.F90')

! ******************************************************************************
! Region-dependent values
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pRegion    => regions(iReg)
    pMixtInput => pRegion%mixtInput
    pSpecInput => pRegion%specInput

! ==============================================================================
!   Gas model
! ==============================================================================

! ------------------------------------------------------------------------------
!   Mixture of thermally and perfect gases: Must not have any discrete species
! ------------------------------------------------------------------------------

    IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_TCPERF ) THEN 
      IF ( global%specUsed .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Species module must be active.')
      ELSE 
#ifdef SPEC
        IF ( pSpecInput%nSpecies <= 1 ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                         'Must have at least one species.')        
        END IF ! pSpecInput%nSpecies
        
        DO iSpec = 1,pSpecInput%nSpecies
          pSpecType => pRegion%specInput%specType(iSpec)
          
          IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Can only have gaseous species.')
          END IF ! pSpecType%discreteFlag
        END DO ! iSpec
#else
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Can only be used with species module.')
#endif      
      END IF ! global%specUsed
    
! ------------------------------------------------------------------------------
!   Pseudogas: Must not have only discrete species
! ------------------------------------------------------------------------------

    ELSE IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_PSEUDO ) THEN 
      IF ( global%specUsed .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Species module must be active.')
      ELSE 
#ifdef SPEC
        IF ( pSpecInput%nSpecies <= 1 ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                         'Must have at least one species.')        
        END IF ! pSpecInput%nSpecies
        
        contCntr = 0
        
        DO iSpec = 1,pSpecInput%nSpecies
          pSpecType => pRegion%specInput%specType(iSpec)
                                                            
          IF ( pSpecType%discreteFlag .EQV. .FALSE. ) THEN 
            contCntr = contCntr + 1
          END IF ! pSpecType%discreteFlag
        END DO ! iSpec
        
        IF ( contCntr == 0 ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                         'Must have at least one gaseous species.')        
        END IF ! contCntr        
#else
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Can only be used with species module.')
#endif      
      END IF ! global%specUsed

! ------------------------------------------------------------------------------
!   Gas-liquid mixture: Must not have any discrete species
! ------------------------------------------------------------------------------

    ELSE IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ ) THEN 
      IF ( global%specUsed .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Species module must be active.')
      ELSE 
#ifdef SPEC
        IF ( pSpecInput%nSpecies <= 1 ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                         'Must have two species.')        
        END IF ! pSpecInput%nSpecies
        
        DO iSpec = 1,pSpecInput%nSpecies
          pSpecType => pRegion%specInput%specType(iSpec)
          
          IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                           'Can only have gaseous species.')
          END IF ! pSpecType%discreteFlag
        END DO ! iSpec   
#else
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, & 
                       'Can only be used with species module.')
#endif      
      END IF ! global%specUsed
    END IF ! pMixtInput%gasModel

! ==============================================================================
!   Mass fluxes
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      pMixtInput%indMfMixt = 1
    ELSE
      pMixtInput%indMfMixt = 0
    END IF ! global%specUsed

! ==============================================================================
!   Substantial derivative
! ==============================================================================
    
#ifdef SPEC
    IF ( global%specUsed .EQV. .TRUE. ) THEN
      DO iSpec = 1,pSpecInput%nSpecies
        IF ( pSpecInput%specType(iSpec)%velocityMethod == &
             SPEC_METHV_EQEUL ) THEN
          pMixtInput%indSd = 1
        END IF ! velocityMethod
      END DO ! iSpec
    END IF ! global%specUsed
#endif

! ==============================================================================
!   Source terms
! ==============================================================================

#ifdef SPEC
    IF ( global%specUsed .EQV. .TRUE. ) THEN
      IF ( (pMixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ) .AND. & 
           (pSpecInput%sourceFlag .EQV. .TRUE.) ) THEN         
        DO iSpec = 1,pSpecInput%nSpecies
          IF ( pSpecInput%specType(iSpec)%sourceType /= &
               SPEC_SOURCE_TYPE_CAVI ) THEN
            CALL ErrorStop(global,ERR_SPEC_SOURCE_TYPE_INVALID,__LINE__)            
          END IF ! sourceType
        END DO ! iSpec
      END IF ! pMixtInput%gasModel
    END IF ! global%specUsed    
#endif
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetDerivedUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetDerivedUserInput.F90,v $
! Revision 1.12  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.9  2006/03/30 20:48:12  haselbac
! Added check of src term for gasliq model
!
! Revision 1.8  2006/03/26 20:21:42  haselbac
! Added support for GL model
!
! Revision 1.7  2005/11/14 16:55:41  haselbac
! Added support for pseudo-gas model
!
! Revision 1.6  2005/11/10 02:55:56  haselbac
! Added check for gas model
!
! Revision 1.5  2005/03/31 16:52:04  haselbac
! Removed setting of nSd
!
! Revision 1.4  2004/10/19 19:25:30  haselbac
! Bug fix: Should only access species data if SPEC == 1
!
! Revision 1.3  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.2  2004/07/28 15:29:19  jferry
! created global variable for spec use
!
! Revision 1.1  2004/01/29 22:56:19  haselbac
! Initial revision
!
! ******************************************************************************







