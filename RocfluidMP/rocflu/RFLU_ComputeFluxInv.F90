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
! Purpose: Wrapper function for computing inviscid fluxes of mixture.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region data
!   fluxPart	Part of flux (central or dissipation or both)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ComputeFluxInv.F90,v 1.3 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************


SUBROUTINE RFLU_ComputeFluxInv(pRegion,fluxPart)

  USE ModDataTypes
  USE ModParameters  
  USE ModError  
  USE ModGlobal, ONLY: t_global  
  USE ModDataStruct, ONLY: t_region
  
  USE RFLU_ModAUSMFlux
  USE RFLU_ModHLLCFlux
  USE RFLU_ModRoeFlux

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: fluxPart
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeFluxInv.F90,v $ $Revision: 1.3 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeFluxInv',&
  'RFLU_ComputeFluxInv.F90')
          
! ******************************************************************************
! Compute fluxes
! ******************************************************************************
         
  SELECT CASE ( pRegion%mixtInput%spaceDiscr )
    CASE ( DISCR_UPW_ROE )        
      CALL RFLU_ROE_ComputeFlux(pRegion,fluxPart)
    CASE ( DISCR_UPW_HLLC )
      IF ( fluxPart == FLUX_PART_CENTRAL .OR. fluxPart == FLUX_PART_BOTH ) THEN
        CALL RFLU_HLLC_ComputeFlux(pRegion)      
      END IF ! fluxPart
    CASE ( DISCR_UPW_AUSMPLUS )     
      IF ( fluxPart == FLUX_PART_CENTRAL .OR. fluxPart == FLUX_PART_BOTH ) THEN     
        CALL RFLU_AUSM_ComputeFlux(pRegion)
      END IF ! fluxPart
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END SELECT ! pMixtInput%spaceDiscr

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeFluxInv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeFluxInv.F90,v $
! Revision 1.3  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/05/01 21:02:58  haselbac
! Initial revision
!
! ******************************************************************************







