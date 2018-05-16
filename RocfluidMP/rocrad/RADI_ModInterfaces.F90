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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_ModInterfaces.F90,v 1.6 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE RADI_ModInterfaces

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Interfaces to external code
! =============================================================================

  SUBROUTINE RADI_CalcEffTemp( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_CalcEffTemp

  SUBROUTINE RADI_CheckParamInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_CheckParamInput

  SUBROUTINE RADI_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_DerivedInputValues

  SUBROUTINE RADI_DiffRadFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_DiffRadFlux

  SUBROUTINE RADI_DiffRadFluxPatch( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_DiffRadFluxPatch

  SUBROUTINE RADI_DiffRadIntens( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_DiffRadIntens

  SUBROUTINE RADI_ExtinctionCoef( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_ExtinctionCoef

  SUBROUTINE RADI_FlimDiffFlux( region )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region) :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
  END SUBROUTINE RADI_FlimDiffFlux

  SUBROUTINE RADI_FlimDiffFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FlimDiffFluxPatch

  SUBROUTINE RADI_FlimEmsInit( region,iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE RADI_FlimEmsInit

  SUBROUTINE RADI_FlimEmsUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimEmsUpdate

  SUBROUTINE RADI_FlimRkInit( region, iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE RADI_FlimRkInit

  SUBROUTINE RADI_FlimRkUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimRkUpdate

  SUBROUTINE RADI_FlimSourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimSourceTerms

  SUBROUTINE RADI_FluxLimiter( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FluxLimiter

  SUBROUTINE RADI_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_InitInputValues

  SUBROUTINE RADI_MixtSourceTermsFlim( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_MixtSourceTermsFlim

  SUBROUTINE RADI_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_ReadInputFile

  SUBROUTINE RADI_ReadRadiSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_ReadRadiSection

  SUBROUTINE RADI_PeulSourceTermsFlim( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_PeulSourceTermsFlim

  SUBROUTINE RADI_PlagSourceTermsFlim( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_PlagSourceTermsFlim

#ifdef RFLO
! =============================================================================
! Rocflo-specific routines
! =============================================================================

  SUBROUTINE RADI_FloFlimBcondDiffuse( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimBcondDiffuse

  SUBROUTINE RADI_FloFlimBcondInjection( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimBcondInjection

  SUBROUTINE RADI_FloFlimBcondSymmetry( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimBcondSymmetry 

  SUBROUTINE RADI_FloFlimBcondZeroGrad( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimBcondZeroGrad

  SUBROUTINE RADI_FloFlimCentFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FloFlimCentFlux

  SUBROUTINE RADI_FloFlimCentFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimCentFluxPatch

  SUBROUTINE RADI_FloFlimCentralDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FloFlimCentralDissipation

  SUBROUTINE RADI_FloFlimCorrCornEdgeCells( region,patch,bcType )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
    INTEGER        :: bcType
  END SUBROUTINE RADI_FloFlimCorrCornEdgeCells

  SUBROUTINE RADI_FloFlimExchangeDummyConf( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RADI_FloFlimExchangeDummyConf

  SUBROUTINE RADI_FloFlimExchCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_FloFlimExchCornEdgeCells

  SUBROUTINE RADI_FloFlimRecvCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_FloFlimRecvCornEdgeCells

  SUBROUTINE RADI_FloFlimRecvDummyVals( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RADI_FloFlimRecvDummyVals

  SUBROUTINE RADI_FloFlimSendCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_FloFlimSendCornEdgeCells

  SUBROUTINE RADI_FloFlimSendDummyConf( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RADI_FloFlimSendDummyConf

  SUBROUTINE RADI_FloFlimSetCornEdgeCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FloFlimSetCornEdgeCells
#endif

  END INTERFACE

END MODULE RADI_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_ModInterfaces.F90,v $
! Revision 1.6  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/09/30 17:11:30  wasistho
! prepared for full FLD radiation model
!
! Revision 1.3  2004/09/18 17:42:11  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.2  2003/07/17 01:12:11  wasistho
! initial activation rocrad
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************






