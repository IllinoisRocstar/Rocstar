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
! $Id: PEUL_ModInterfaces.F90,v 1.10 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE PEUL_ModInterfaces

  IMPLICIT NONE

  INTERFACE

! These are PUBLIC (i.e., in ModInterfacesEulerian)
!
! PEUL_AllocateDataBuffers, PEUL_AllocateMemory,
! PEUL_boundaryConditionsRecv, PEUL_boundaryConditionsSend,
! PEUL_BoundaryConditionsSet,
! PEUL_CentralDissipation, PEUL_ClearSendRequests, PEUL_ConvectiveFluxes,
! PEUL_EnforcePositivity,
! PEUL_InitSolution, PEUL_PrintUserInput, PEUL_ReadSolution,
! PEUL_ResidualSmoothing, PEUL_ResidualSmoothingCoeffs,
! PEUL_SourceTerms,
! PEUL_SpectralRadii PEUL_UserInput, PEUL_WriteSolution

! These are PRIVATE (i.e., not in ModInterfacesEulerian)

  SUBROUTINE PEUL_BcondInflow( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_BcondInflow

  SUBROUTINE PEUL_BcondInjection( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_BcondInjection

  SUBROUTINE PEUL_BcondOutflow( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_BcondOutflow

  SUBROUTINE PEUL_BcondSlipWall( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_BcondSlipWall

  SUBROUTINE PEUL_BcondSymmetry( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_BcondSymmetry

  SUBROUTINE PEUL_CentralFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_CentralFlux

  SUBROUTINE PEUL_CentralFluxPatch( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
  END SUBROUTINE PEUL_CentralFluxPatch

  SUBROUTINE PEUL_CorrectCornerEdgeCells( region,patch,bcType )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_patch),  INTENT(IN)    :: patch
    INTEGER,        INTENT(IN)    :: bcType
  END SUBROUTINE PEUL_CorrectCornerEdgeCells

  SUBROUTINE PEUL_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_DerivedInputValues

  SUBROUTINE PEUL_ExchangeCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_ExchangeCornerEdgeCells

  SUBROUTINE PEUL_ExchangeDummyConf( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_region), INTENT(IN)    :: regionSrc
    TYPE(t_patch),  INTENT(INOUT) :: patch
    TYPE(t_patch),  INTENT(IN)    :: patchSrc
  END SUBROUTINE PEUL_ExchangeDummyConf

  SUBROUTINE PEUL_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_InitInputValues

  SUBROUTINE PEUL_ReadBcFarfSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadBcFarfSection

  SUBROUTINE PEUL_ReadBcInflowSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadBcInflowSection

  SUBROUTINE PEUL_ReadBcInjectSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadBcInjectSection

  SUBROUTINE PEUL_ReadBcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadBcInputFile

  SUBROUTINE PEUL_ReadConPartPtypeSection( regions,brbeg,brend,iPtype )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: brbeg,brend,iPtype
  END SUBROUTINE PEUL_ReadConPartPtypeSection

  SUBROUTINE PEUL_ReadConPartSection( regions,nPtypes,brbeg,brend )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: nPtypes
    INTEGER, INTENT(OUT)    :: brbeg,brend
  END SUBROUTINE PEUL_ReadConPartSection

  SUBROUTINE PEUL_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadInputFile

  SUBROUTINE PEUL_ReceiveCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_ReceiveCornerEdgeCells

  SUBROUTINE PEUL_ReceiveDummyVals( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    TYPE(t_region), INTENT(IN)    :: regionSrc
    TYPE(t_patch),  INTENT(INOUT) :: patch
    TYPE(t_patch),  INTENT(IN)    :: patchSrc
  END SUBROUTINE PEUL_ReceiveDummyVals

  SUBROUTINE PEUL_SendCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_SendCornerEdgeCells

  SUBROUTINE PEUL_SendDummyConf( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region, regionSrc
    TYPE(t_patch),  INTENT(INOUT) :: patch
  END SUBROUTINE PEUL_SendDummyConf

  SUBROUTINE PEUL_SetCornerEdgeCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_SetCornerEdgeCells

  SUBROUTINE PEUL_SourceEqEul( region,ipt )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: ipt
  END SUBROUTINE PEUL_SourceEqEul

  END INTERFACE

END MODULE PEUL_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ModInterfaces.F90,v $
! Revision 1.10  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.7  2004/03/02 21:42:47  jferry
! Added clipping options and corner and edge cell updates
!
! Revision 1.6  2003/04/09 15:10:22  jferry
! added slip wall boundary condition
!
! Revision 1.5  2003/04/09 14:30:30  fnajjar
! Added Interfaces for Multi-region and MPI-based routines
!
! Revision 1.4  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.3  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
! Revision 1.2  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.1.1.1  2002/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************






