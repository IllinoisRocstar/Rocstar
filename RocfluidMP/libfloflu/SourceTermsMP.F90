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
! Purpose: add source terms to the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: pRegion%levels%mixt%rhs = complete right-hand side (residual).
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: SourceTermsMP.F90,v 1.8 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SourceTermsMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : SourceTerms
  USE ModError
  USE ModParameters

#ifdef RFLU
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons
#endif

#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_SourceTerms
#endif
#ifdef RFLU
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_CorrectMixtProperties
#endif
#ifdef SPEC
  USE SPEC_RFLU_ModChemistry, ONLY: SPEC_RFLU_IntegrateChemSrcTerm
#endif
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_SourceTerms
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_SourceTerms
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_SourceTerms
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_RansSourceTerms
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_RFLU_SourceTerms_GL 
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'SourceTermsMP',&
  'SourceTermsMP.F90' )

! add source terms
  
  CALL SourceTerms( region )

#ifdef RFLU
  pRegion => region
  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
  
#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN 
    CALL PLAG_RFLU_CorrectMixtProperties(pRegion)
  END IF ! global%plagUsed
#endif  
#endif

#ifdef INRT
  IF (global%inrtUsed) THEN
    CALL INRT_SourceTerms( region )
  ENDIF
#endif
#ifdef PEUL
  IF (global%peulUsed) THEN
    CALL PEUL_SourceTerms( region )
  ENDIF
#endif
#ifdef RADI
  CALL RADI_SourceTerms( region )
#endif
#ifdef PERI
  CALL PERI_SourceTerms( region )
#endif
#ifdef TURB
  IF ((region%mixtInput%flowModel == FLOW_NAVST) .AND. &
      (region%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
    CALL TURB_RansSourceTerms( region )
  ENDIF
#endif


#ifdef RFLU
  pRegion => region

  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    IF ( pRegion%specInput%sourceFlag .EQV. .TRUE. ) THEN   
   
! ------------------------------------------------------------------------------
!     Cavitation source term for gas-liquid-vapor mixture model
! ------------------------------------------------------------------------------   
   
      IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ ) THEN
        CALL SPEC_RFLU_SourceTerms_GL(pRegion)
      END IF ! pRegion%mixtInput%gasModel

! ------------------------------------------------------------------------------
!     Combustion source term for burning-crack simulations
! ------------------------------------------------------------------------------   

! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that no
!             longer call source term integration function directly in 
!             RFLU_TimeStepping.F90.
!      CALL SPEC_RFLU_IntegrateChemSrcTerm(pRegion,0)
! END TEMPORARY
    END IF ! pRegion%specInput%sourceFlag
  END IF ! global%specUsed
#endif
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE SourceTermsMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SourceTermsMP.F90,v $
! Revision 1.8  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/03/30 20:47:22  haselbac
! Clean-up of source terms for species
!
! Revision 1.5  2006/03/26 20:21:22  haselbac
! Added call to GL src term routine
!
! Revision 1.4  2005/11/30 22:05:25  fnajjar
! Added call to PLAG_RFLU_CorrectMixtProperties
!
! Revision 1.3  2005/10/05 13:48:57  haselbac
! Bug fix: Enclosed call to chem src term within IF
!
! Revision 1.2  2005/06/06 14:23:03  haselbac
! Adapted to Lucas changes
!
! Revision 1.1  2004/12/01 16:51:28  haselbac
! Initial revision after changing case
!
! Revision 1.15  2004/05/03 15:09:41  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.14  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.13  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.12  2004/01/31 03:56:19  haselbac
! Added RFLU state vector conversion routines
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/10/03 20:13:02  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.7  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.6  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.5  2003/08/28 20:33:00  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.4  2003/07/17 00:55:20  wasistho
! initial activation rocrad
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/05 01:59:33  wasistho
! install ROCPERI
!
! Revision 1.1  2003/03/28 19:47:43  fnajjar
! Initial import for RocfluidMP
!
! ******************************************************************************







