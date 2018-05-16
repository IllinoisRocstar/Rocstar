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
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesTurbulence.F90,v 1.23 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesTurbulence

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! turbulence
! =============================================================================

  SUBROUTINE TURB_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE TURB_AllocateMemory

  SUBROUTINE TURB_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE TURB_BuildVersionString

  SUBROUTINE TURB_CalcMetrics( regions, isInit )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER        :: isInit
  END SUBROUTINE TURB_CalcMetrics

  SUBROUTINE TURB_EmsInit( region,iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE TURB_EmsInit

  SUBROUTINE TURB_InitSolution( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_InitSolution

  SUBROUTINE TURB_LesCommunication( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_LesCommunication

  SUBROUTINE TURB_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_PrintUserInput

  SUBROUTINE TURB_RansConvectiveFluxes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansConvectiveFluxes

  SUBROUTINE TURB_RansNumericalDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansNumericalDissipation

  SUBROUTINE TURB_RansSourceTerms( region ) 
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansSourceTerms

  SUBROUTINE TURB_RansZeroDummyCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansZeroDummyCells

  SUBROUTINE TURB_RkInit( region, iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE TURB_RkInit

  SUBROUTINE TURB_SolutionUpdate( region,iStage,ibc,iec )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage, ibc, iec
  END SUBROUTINE TURB_SolutionUpdate

  SUBROUTINE TURB_StatMapping( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE TURB_StatMapping

  SUBROUTINE TURB_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_UserInput

  SUBROUTINE TURB_CoViscousFluxes( region )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
  END SUBROUTINE TURB_CoViscousFluxes

#ifdef RFLO
! =============================================================================
! Rocflo-specific routines
! =============================================================================

  SUBROUTINE TURB_RFLO_RansAllocDataBuffers( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_RFLO_RansAllocDataBuffers

  SUBROUTINE TURB_RFLO_RansBndConditionsRecv( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_RFLO_RansBndConditionsRecv

  SUBROUTINE TURB_RFLO_RansBndConditionsSend( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_RFLO_RansBndConditionsSend

  SUBROUTINE TURB_RFLO_RansBndConditionsSet( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_RFLO_RansBndConditionsSet

  SUBROUTINE TURB_RFLO_RansClearSendRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_RFLO_RansClearSendRequests

  SUBROUTINE TURB_RFLO_RansResSmoothing( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RFLO_RansResSmoothing

  SUBROUTINE TURB_RFLO_RansResSmoothingCoeff( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RFLO_RansResSmoothingCoeff

  SUBROUTINE TURB_RFLO_RansSpectralRadii( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RFLO_RansSpectralRadii

  SUBROUTINE TURB_RFLO_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_RFLO_ReadSolution

  SUBROUTINE TURB_RFLO_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_RFLO_ReadSolutionRegion

  SUBROUTINE TURB_RFLO_WriteSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_RFLO_WriteSolution
#endif

#ifdef RFLU
! =============================================================================
! Rocflu-specific routines
! =============================================================================

  SUBROUTINE TURB_RFLU_DeallocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_RFLU_DeallocateMemory

  SUBROUTINE TURB_RFLU_ReadSolutionASCII( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_RFLU_ReadSolutionASCII

  SUBROUTINE TURB_RFLU_ReadSolutionBinary( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_RFLU_ReadSolutionBinary

  SUBROUTINE TURB_RFLU_WriteSolutionASCII( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_RFLU_WriteSolutionASCII

  SUBROUTINE TURB_RFLU_WriteSolutionBinary( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_RFLU_WriteSolutionBinary
#endif

  END INTERFACE

END MODULE ModInterfacesTurbulence

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesTurbulence.F90,v $
! Revision 1.23  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.22  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.21  2004/12/01 00:08:40  wasistho
! added BuildVersionString
!
! Revision 1.20  2004/11/17 23:43:50  wasistho
! used generic RK-update for rocturb
!
! Revision 1.19  2004/07/03 01:11:02  wasistho
! set it back to the previous version
!
! Revision 1.18  2004/07/03 00:51:11  wasistho
! make argument types the same for TURB_coViscousFluxes btw FLO and FLU
!
! Revision 1.17  2004/06/19 03:29:33  wasistho
! removed argument iReg in TURB_InitSolution
!
! Revision 1.16  2004/03/27 02:20:14  wasistho
! compiled with Rocflu
!
! Revision 1.15  2004/03/20 00:26:52  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.14  2004/03/19 02:40:40  wasistho
! renamed TURB_RFLO_RansZeroDummyCells to TURB_RansZeroDummyCells
!
! Revision 1.13  2004/03/13 03:07:43  wasistho
! get rid of flo/flu identifier in TURB_Co.. routine names
!
! Revision 1.12  2004/03/11 03:30:09  wasistho
! changed rocturb nomenclature
!
! Revision 1.11  2004/03/08 23:33:08  wasistho
! changed turb file/routine name
!
! Revision 1.10  2004/02/26 21:18:23  wasistho
! added TURB_rkInit, TURB_lesComm, TURB_emsInit
!
! Revision 1.9  2004/02/11 03:23:17  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.8  2004/02/07 00:55:52  wasistho
! added TURB_ReadSolutionRegion in interface module
!
! Revision 1.6  2003/10/27 04:50:09  wasistho
! added RaNS upwind schemes
!
! Revision 1.5  2003/10/16 20:19:28  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
! Revision 1.4  2003/10/03 20:14:21  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.3  2003/08/06 15:55:08  wasistho
! added CalcVortic and SolutionUpdate for vorticities
!
! Revision 1.2  2003/07/22 02:53:50  wasistho
! prepare accurate rocturb restart
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






