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
! $Id: RFLO_ModInterfacesSolver.F90,v 1.30 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE RFLO_ModInterfacesSolver

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! ROCFLO solver
! =============================================================================

  SUBROUTINE RFLO_AllocateDataBuffers( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RFLO_AllocateDataBuffers

  SUBROUTINE RFLO_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_AllocateMemory

  SUBROUTINE RFLO_C2eAvgCoeffs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_C2eAvgCoeffs

  SUBROUTINE RFLO_C2eAvgCoeffsDegec( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_C2eAvgCoeffsDegec

  SUBROUTINE RFLO_C2fAvgCoeffs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_C2fAvgCoeffs

  SUBROUTINE RFLO_C2fAvgCoeffsDegec( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_C2fAvgCoeffsDegec

  SUBROUTINE RFLO_C2fAvgCoeffsDummy( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_C2fAvgCoeffsDummy

  SUBROUTINE RFLO_C2fAvgCoeffsDummyConn( region,lbound,idir,jdir,kdir, &
                             indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: lbound,idir,jdir,kdir
    INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  END SUBROUTINE RFLO_C2fAvgCoeffsDummyConn

  SUBROUTINE RFLO_C2fAvgCoeffsDummyPhys( region,lbound,idir,jdir,kdir, &
                             indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: lbound,idir,jdir,kdir
    INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  END SUBROUTINE RFLO_C2fAvgCoeffsDummyPhys

  SUBROUTINE RFLO_C2fAvgCoeffsPatch( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_C2fAvgCoeffsPatch

  SUBROUTINE RFLO_CalcControlVolumes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcControlVolumes

  SUBROUTINE RFLO_CalcFaceVectors( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcFaceVectors

  SUBROUTINE RFLO_CalcForces( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcForces

  SUBROUTINE RFLO_CalcGridSpeeds( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcGridSpeeds

  SUBROUTINE RFLO_CalcMassFlow( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcMassFlow

  SUBROUTINE RFLO_CalcThrust( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcThrust

  SUBROUTINE RFLO_CalcTotalMass( region,mass )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL) :: mass
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcTotalMass

  SUBROUTINE RFLO_CentralFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CentralFlux

  SUBROUTINE RFLO_CentralFluxPatch( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_CentralFluxPatch

  SUBROUTINE RFLO_CentralDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CentralDissipation

  SUBROUTINE RFLO_CheckBcInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CheckBcInput

  SUBROUTINE RFLO_CheckDerivedUserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CheckDerivedUserInput

  SUBROUTINE RFLO_CheckRegionFaces( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CheckRegionFaces

  SUBROUTINE RFLO_CheckMetrics( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CheckMetrics

  SUBROUTINE RFLO_CheckMinimumCells( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CheckMinimumCells

  SUBROUTINE RFLO_CheckUserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CheckUserInput

  SUBROUTINE RFLO_ClearSendRequests( regions,iReg,geometry )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
    LOGICAL :: geometry
  END SUBROUTINE RFLO_ClearSendRequests

#ifndef GENX
  SUBROUTINE RFLO_ComputeIntegralValues(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ComputeIntegralValues
#else
  SUBROUTINE RFLO_ComputeIntegralValues(regions,integ)
    USE ModRocstar ! To access MAN_INTEG_SIZE
    USE ModDataStruct, ONLY: t_region
    DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ComputeIntegralValues
#endif

  SUBROUTINE RFLO_CopyBoundaryData( global,patchPrev,patch )
    USE ModGlobal, ONLY     : t_global
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER  :: patchPrev, patch
  END SUBROUTINE RFLO_CopyBoundaryData

  SUBROUTINE RFLO_CopyGeometryDummy( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CopyGeometryDummy

  SUBROUTINE RFLO_CopyTopologyLevels( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_CopyTopologyLevels

  SUBROUTINE RFLO_CopyVectorCorners( iLev,region,vec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: vec(:)
  END SUBROUTINE RFLO_CopyVectorCorners

  SUBROUTINE RFLO_CopyVectorEdges( iLev,region,vec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: vec(:)
  END SUBROUTINE RFLO_CopyVectorEdges

  SUBROUTINE RFLO_CopyVectorPatches( iLev,region,vec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: vec(:)
  END SUBROUTINE RFLO_CopyVectorPatches

  SUBROUTINE RFLO_CopyMatrixCorners( iLev,region,mat )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: mat(:,:)
  END SUBROUTINE RFLO_CopyMatrixCorners

  SUBROUTINE RFLO_CopyMatrixEdges( iLev,region,mat )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: mat(:,:)
  END SUBROUTINE RFLO_CopyMatrixEdges

  SUBROUTINE RFLO_CopyMatrixPatches( iLev,region,mat )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iLev
    REAL(RFREAL), POINTER :: mat(:,:)
  END SUBROUTINE RFLO_CopyMatrixPatches

  SUBROUTINE RFLO_CorrectCornerEdgeCells( region,patch,bcType )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: bcType
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_CorrectCornerEdgeCells

  SUBROUTINE RFLO_DoMemoryAllocation( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE (t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DoMemoryAllocation

  SUBROUTINE RFLO_DualMultigrid( dTimeSystem,regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL) :: dTimeSystem
    TYPE (t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DualMultigrid

  SUBROUTINE RFLO_DualTimeStepping( dTimeSystem,regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL) :: dTimeSystem
    TYPE (t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DualTimeStepping

  SUBROUTINE RFLO_DualTstInit( regions,timeLevel )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: timeLevel
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DualTstInit

  SUBROUTINE RFLO_DualTstPredict( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_DualTstPredict

  SUBROUTINE RFLO_DualTstSterm( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_DualTstSterm

  SUBROUTINE RFLO_DualTstShift( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_DualTstShift

  SUBROUTINE RFLO_EndFlowSolver( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_EndFlowSolver

  SUBROUTINE RFLO_ExchangeCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RFLO_ExchangeCornerEdgeCells

  SUBROUTINE RFLO_ExchangeDummyConf( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeDummyConf

  SUBROUTINE RFLO_ExchangeDummyInt( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeDummyInt

  SUBROUTINE RFLO_ExchangeDummyIreg( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeDummyIreg

  SUBROUTINE RFLO_ExchangeGeometry( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ExchangeGeometry

  SUBROUTINE RFLO_ExchangeGeometryCopy( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeGeometryCopy

  SUBROUTINE RFLO_ExchangeGeometryLevels( region,iPatch )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iPatch
  END SUBROUTINE RFLO_ExchangeGeometryLevels

  SUBROUTINE RFLO_ExchangeGeometryPrepare( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ExchangeGeometryPrepare

  SUBROUTINE RFLO_ExchangeGeometryRecv( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeGeometryRecv

  SUBROUTINE RFLO_ExchangeGeometrySend( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_ExchangeGeometrySend

  SUBROUTINE RFLO_ExtrapolateGeometry( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ExtrapolateGeometry

  SUBROUTINE RFLO_FindSourceCell( regions,iReg,iLev,ic,jc,kc,icell, &
                                  found,rotate,iRegSrc )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg, iLev, ic, jc, kc, icell, iRegSrc
    LOGICAL :: found, rotate
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_FindSourceCell

  SUBROUTINE RFLO_FindSourceCellInvert( regions,iReg,iLev,ic,jc,kc, &
                                        icell,found,rotate,iRegSrc )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg, iLev, ic, jc, kc, icell, iRegSrc
    LOGICAL :: found, rotate
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_FindSourceCellInvert

  SUBROUTINE RFLO_SourceCell( region,regionSrc,patch,patchSrc, &
                              iLev, ic,jc,kc,icell,found )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iLev, ic, jc, kc, icell
    LOGICAL :: found
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch), POINTER :: patch, patchSrc
  END SUBROUTINE RFLO_SourceCell

  SUBROUTINE RFLO_FindSourceRegions( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_FindSourceRegions

  SUBROUTINE RFLO_FindSourcePatches( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_FindSourcePatches

  SUBROUTINE RFLO_FindThrustPatches( region,iReg )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_FindThrustPatches

  SUBROUTINE RFLO_InitAvgCoeffs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_InitAvgCoeffs

#ifdef GENX
  SUBROUTINE RFLO_InitFlowSolver( globalGenx,initialTime,communicator, &
                                  genxHandle,inSurf,inVol,obtain_attribute )
    USE ModRocstar, ONLY : t_globalGenx
    CHARACTER(*), INTENT(in) :: inSurf, inVol
    DOUBLE PRECISION, INTENT(in) :: initialTime
    INTEGER, INTENT(in)  :: communicator, genxHandle, obtain_attribute
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLO_InitFlowSolver

  SUBROUTINE RFLO_FlowSolver( globalGenx,timeSystem,dTimeSystem,genxHandleBc, &
                              genxHandleGm )
    USE ModRocstar, ONLY : t_globalGenx
    INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm
    DOUBLE PRECISION, INTENT(in) :: timeSystem, dTimeSystem
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLO_FlowSolver
#else
  SUBROUTINE RFLO_InitFlowSolver( casename,verbLevel,global,regions )
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal,     ONLY : t_global
    CHARACTER(*) :: casename
    INTEGER      :: verbLevel
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_InitFlowSolver

  SUBROUTINE RFLO_FlowSolver( dTimeSystem,dIterSystem,regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL) :: dTimeSystem
    INTEGER      :: dIterSystem
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_FlowSolver
#endif

  SUBROUTINE RFLO_GetFlowSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_GetFlowSolution

  SUBROUTINE RFLO_GetGeometry( regions,iread )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iread
  END SUBROUTINE RFLO_GetGeometry

  SUBROUTINE RFLO_GetUserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_GetUserInput

  SUBROUTINE RFLO_InterpolToFinerLevel( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_InterpolToFinerLevel

  SUBROUTINE RFLO_InitGridProcedures( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_InitGridProcedures

  SUBROUTINE RFLO_LimiterReference( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_LimiterReference

  SUBROUTINE RFLO_MapRegionsProcessors( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MapRegionsProcessors

  SUBROUTINE RFLO_MinimumTimeStep( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MinimumTimeStep

  SUBROUTINE RFLO_MirrorGeometry( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_MirrorGeometry

  SUBROUTINE RFLO_MoveGridBlocks( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MoveGridBlocks

  SUBROUTINE RFLO_MoveGridGlobal( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MoveGridGlobal

  SUBROUTINE RFLO_MoveGridInterfaces( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MoveGridInterfaces

  SUBROUTINE RFLO_MoveGridSurfaces( regions,someMoved )
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: someMoved
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_MoveGridSurfaces

  SUBROUTINE RFLO_Multigrid( dIterSystem,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: dIterSystem
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_Multigrid

  SUBROUTINE RFLO_NewGrid( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_NewGrid

  SUBROUTINE RFLO_OpenConverFile( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_OpenConverFile

  SUBROUTINE RFLO_OpenProbeFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_OpenProbeFile

  SUBROUTINE RFLO_OpenThrustFile( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_OpenThrustFile

  SUBROUTINE RFLO_PrintUserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_PrintUserInput

  SUBROUTINE RFLO_ReadBcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcInputFile

  SUBROUTINE RFLO_ReadBcFromFile( global,fname,patch )
    USE ModGlobal, ONLY   : t_global
    USE ModBndPatch, ONLY : t_patch
    CHARACTER(*) :: fname
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER  :: patch
  END SUBROUTINE RFLO_ReadBcFromFile

  SUBROUTINE RFLO_ReadBcFarfSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcFarfSection

  SUBROUTINE RFLO_ReadBcNoslipSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcNoslipSection

  SUBROUTINE RFLO_ReadBcInflowTotAngSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcInflowTotAngSection

  SUBROUTINE RFLO_ReadBcInflowVelSection( regions,bcTitle )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: bcTitle
  END SUBROUTINE RFLO_ReadBcInflowVelSection

  SUBROUTINE RFLO_ReadBcInjectMrateSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcInjectMrateSection

  SUBROUTINE RFLO_ReadBcInjectAPNSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcInjectAPNSection

  SUBROUTINE RFLO_ReadBcOutflowSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcOutflowSection

  SUBROUTINE RFLO_ReadBcSlipWallSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadBcSlipWallSection

  SUBROUTINE RFLO_ReadRegionMapSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_ReadRegionMapSection

  SUBROUTINE RFLO_ReadTbcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadTbcInputFile

  SUBROUTINE RFLO_ReadTbcSection( regions,tbcType )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: tbcType
  END SUBROUTINE RFLO_ReadTbcSection

  SUBROUTINE RFLO_ReceiveCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RFLO_ReceiveCornerEdgeCells

  SUBROUTINE RFLO_ReceiveDummyVals( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ReceiveDummyVals

  SUBROUTINE RFLO_ResidualNorm( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ResidualNorm

  SUBROUTINE RFLO_ResidualSmoothing( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ResidualSmoothing

  SUBROUTINE RFLO_ResidualSmoothingCoeffs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ResidualSmoothingCoeffs

  SUBROUTINE RFLO_RoeDissipFirst( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_RoeDissipFirst

  SUBROUTINE RFLO_RoeDissipSecond( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_RoeDissipSecond

  SUBROUTINE RFLO_RoeFluxFirst( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_RoeFluxFirst

  SUBROUTINE RFLO_RoeFluxSecond( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_RoeFluxSecond

  SUBROUTINE RFLO_RoeFluxPatch( region,patch )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_RoeFluxPatch

  SUBROUTINE RFLO_SendCornerEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RFLO_SendCornerEdgeCells

  SUBROUTINE RFLO_SendDummyConf( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_SendDummyConf

  SUBROUTINE RFLO_SendDummyInt( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_SendDummyInt

  SUBROUTINE RFLO_SendDummyIreg( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_SendDummyIreg

  SUBROUTINE RFLO_SetCornerEdgeCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_SetCornerEdgeCells

  SUBROUTINE RFLO_TimeStepInviscid( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_TimeStepInviscid

  SUBROUTINE RFLO_TimeStepViscous( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_TimeStepViscous

  SUBROUTINE RFLO_TimeStepping( dTimeSystem,dIterSystem,regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL) :: dTimeSystem
    INTEGER      :: dIterSystem
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_TimeStepping

  SUBROUTINE RFLO_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_UserInput

  SUBROUTINE RFLO_ViscousFlux( region,indxMu,indxTCo,tv )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: indxMu, indxTCo
    REAL(RFREAL), POINTER :: tv(:,:)
  END SUBROUTINE RFLO_ViscousFlux

  SUBROUTINE RFLO_ViscousFluxPatch( region,patch,indxMu,indxTCo,tv )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
    INTEGER        :: indxMu, indxTCo
    REAL(RFREAL), POINTER :: tv(:,:)
  END SUBROUTINE RFLO_ViscousFluxPatch

  END INTERFACE

END MODULE RFLO_ModInterfacesSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModInterfacesSolver.F90,v $
! Revision 1.30  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.29  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.28  2006/03/04 04:31:58  wasistho
! moved RFLO_CalcGridMetrics to rocflo module
!
! Revision 1.27  2006/01/20 06:14:55  wasistho
! added ReadBcInjectMrate and ReadBcInjectAPN
!
! Revision 1.26  2005/11/11 07:16:04  wasistho
! removed RFLO_FindDegeneratCell
!
! Revision 1.25  2005/10/20 06:53:41  wasistho
! removed RFLO_CalcCellCentroids from list
!
! Revision 1.24  2005/05/28 08:04:13  wasistho
! added RFLO_InitGridProcedures
!
! Revision 1.23  2005/05/27 08:09:18  wasistho
! added argument iread in rflo_getgeometry
!
! Revision 1.22  2005/04/28 22:04:50  wasistho
! added RFLO_ReadBcInflow...
!
! Revision 1.21  2005/03/01 16:36:02  wasistho
! added ModGenx in RFLO_ComputeIntegralValues
!
! Revision 1.20  2005/02/26 04:06:25  wasistho
! added RFLO_ComputeIntegralValues
!
! Revision 1.19  2004/12/28 22:50:20  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
! Revision 1.18  2004/12/02 23:28:02  wasistho
! removed entry BuildVersionString
!
! Revision 1.17  2004/11/30 20:10:14  fnajjar
! Included interface for RFLO_CheckDerivedUserInput
!
! Revision 1.16  2004/08/25 07:47:35  wasistho
! added RFLO_C2f/eAvgCoeffsDegec and RFLO_InitAvgCoeffs
!
! Revision 1.15  2004/08/21 00:35:21  wasistho
! added RFLO_findSourceCellInvert and RFLO_findDegeneratCell
!
! Revision 1.14  2004/08/03 22:45:11  wasistho
! added RFLO_c2eAvgCoeffs
!
! Revision 1.13  2004/08/02 23:12:27  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.12  2004/07/30 17:28:27  wasistho
! added routines starting RFLO_c2f...
!
! Revision 1.11  2004/06/29 23:58:08  wasistho
! migrated to Roccom-3
!
! Revision 1.10  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.9  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.8  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.7  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.6  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.5  2003/05/20 20:46:57  jblazek
! Values in edge & corner cells now corrected at noslip and symmetry walls.
!
! Revision 1.4  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.3  2003/02/03 19:20:46  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.2  2003/01/15 22:10:13  jblazek
! Other interfaces to InitFlowSolver and FlowSolver with GENX.
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






