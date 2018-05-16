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
! Purpose: Wrapper for checking validity of variables.
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
! $Id: RFLU_CheckValidityWrapper.F90,v 1.3 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckValidityWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_CheckValidity, &
                           RFLU_CheckValidity_GL

#ifdef PLAG
  USE PLAG_ModCheckVars, ONLY: PLAG_CheckValidity
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

  RCSIdentString = '$RCSfile: RFLU_CheckValidityWrapper.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckValidityWrapper',&
  'RFLU_CheckValidityWrapper.F90')

! ******************************************************************************
! Mixture
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================

    CASE ( FLUID_MODEL_INCOMP )

! ==============================================================================
!   Compressible fluid model
! ==============================================================================

    CASE ( FLUID_MODEL_COMP )
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
        CASE ( GAS_MODEL_TCPERF, & 
               GAS_MODEL_MIXT_TCPERF, & 
               GAS_MODEL_MIXT_TPERF, & 
               GAS_MODEL_MIXT_PSEUDO )
          CALL RFLU_CheckValidity(pRegion)
        CASE ( GAS_MODEL_MIXT_GASLIQ )
          CALL RFLU_CheckValidity_GL(pRegion) 
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

#ifdef PLAG
! ==============================================================================
! Particles
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL PLAG_CheckValidity(pRegion)
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckValidityWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckValidityWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/03/26 20:20:58  haselbac
! Initial revision
!
! ******************************************************************************







