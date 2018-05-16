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
! Purpose: Wrapper for check for posivity of variables
!
! Description: None.
!
! Input:
!   pRegion        Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CheckPositivityWrapper.F90,v 1.6 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckPositivityWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_CheckPositivity, &
                           RFLU_CheckPositivity_GL, &
                           RFLU_CheckValidity_GL 

#ifdef SPEC
  USE ModInterfaces, ONLY: RFLU_ScalarCheckPositivity
#endif

#ifdef PLAG
  USE PLAG_ModCheckVars, ONLY: PLAG_CheckPositivity
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_CheckPositivityWrapper.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckPositivityWrapper',&
  'RFLU_CheckPositivityWrapper.F90')

! ******************************************************************************
! Mixture
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================

     CASE ( FLUID_MODEL_INCOMP )

! ==============================================================================
!   Compressible fluid model. NOTE check state of solution vector.
! ==============================================================================

     CASE ( FLUID_MODEL_COMP )
       SELECT CASE ( pRegion%mixtInput%gasModel ) 
         CASE ( GAS_MODEL_TCPERF, & 
                GAS_MODEL_MIXT_TCPERF, & 
                GAS_MODEL_MIXT_TPERF, & 
                GAS_MODEL_MIXT_PSEUDO )
           CALL RFLU_CheckPositivity(pRegion)
         CASE ( GAS_MODEL_MIXT_GASLIQ )
           CALL RFLU_CheckPositivity_GL(pRegion)
         CASE DEFAULT
           CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
       END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel


! ******************************************************************************
! Multiphysics modules
! ******************************************************************************

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarCheckPositivity(pRegion,FTYPE_SPEC, &
                                    pRegion%specInput%nSpecies, &
                                    pRegion%spec%cv)
  END IF ! global%specUsed
#endif

#ifdef PLAG
! ==============================================================================
! Particles
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL PLAG_CheckPositivity(pRegion)
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckPositivityWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckPositivityWrapper.F90,v $
! Revision 1.6  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/26 20:21:27  haselbac
! Rewrite bcos of GL model, cosmetics
!
! Revision 1.3  2005/12/01 21:54:43  fnajjar
! Added call to PLAG_CheckPositivity
!
! Revision 1.2  2004/07/28 15:29:19  jferry
! created global variable for spec use
!
! Revision 1.1  2004/01/29 22:55:50  haselbac
! Initial revision
!
! ******************************************************************************







