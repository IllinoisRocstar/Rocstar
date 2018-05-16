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
! Purpose: compute source term due to radiation and added to RHS of fluid eqs.
!
! Description: For diffusion approximation, computing source term consists of:
!              - computing radiant energy flux (qr), added to RHS (ROSS/FLDSRC)
!              - computing radiant source term, added to RHS (FLDTRAN)
!              For general method of radiation modeling:
!              - computing radiation intensities I (heart of the method),
!              - computing radiant energy flux (qr) from I and added to RHS.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: rhs of energy equation updated by RADI source term
!
! Notes: effective temperature is computed for all, but for general method
!        and full-FLD it is used merely for comparison and reference, not 
!        employed in the computation.
!
!******************************************************************************
!
! $Id: RADI_SourceTerms.F90,v 1.4 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_SourceTerms( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE RADI_ModInterfaces, ONLY : RADI_DiffRadFlux, RADI_DiffRadIntens
!                                 RADI_ExtinctionCoef, RADI_CalcEffTemp, &
!                                 RADI_FluxLimiter
!                                 RADI_MixtSourceTermsFlim, &
!                                 RADI_GenSolveRTE, RADI_GenRadFlux
  USE ModError
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_SourceTerms',&
  'RADI_SourceTerms.F90' )

! compute radiation quantities and add radiation source term to RHS -----------

#ifdef RFLO
  IF (.NOT. region%mixtInput%radiUsed) GOTO 999

! radiation intensities and radiant energy flux
  IF ((region%radiInput%radiModel == RADI_MODEL_RTEGRAY) .OR. &
      (region%radiInput%radiModel == RADI_MODEL_RTEBAND)) THEN
!    CALL RADI_GenSolveRTE( region )
!    CALL RADI_GenRadFlux( region )
  ELSEIF (region%radiInput%radiModel == RADI_MODEL_ROSS) THEN
!    CALL RADI_ExtinctionCoef( region )
!    CALL RADI_CalcEffTemp( region )
!    CALL RADI_FluxLimiter( region )
    CALL RADI_DiffRadFlux( region )
    CALL RADI_DiffRadIntens( region )
  ELSEIF (region%radiInput%radiModel == RADI_MODEL_FLDSRC) THEN
    CALL RADI_DiffRadFlux( region )
    CALL RADI_DiffRadIntens( region )
  ELSEIF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
!    CALL RADI_MixtSourceTermsFlim( region )
#ifdef PEUL
!    IF (global%peulUsed)
!      CALL RADI_PeulSourceTermsFlim( region )
!    ENDIF
#endif
#ifdef PEUL
!    IF (global%plagUsed)
!      CALL RADI_PlagSourceTermsFlim( region )
!    ENDIF
#endif
  ENDIF ! radiModel
#endif

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_SourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_SourceTerms.F90,v $
! Revision 1.4  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:50  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.4  2004/09/22 01:31:47  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.3  2004/09/18 17:41:35  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.2  2003/07/23 03:14:33  wasistho
! cured baby illness
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************







