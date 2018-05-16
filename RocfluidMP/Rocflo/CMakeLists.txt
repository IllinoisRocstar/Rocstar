cmake_minimum_required (VERSION 2.8)
project(Rocflo Fortran C CXX)
enable_language(Fortran C CXX)
include(FortranCInterface)
#load_cache(../)
set (ROCSTAR TRUE CACHE BOOL "Build as Rocstar Module")
set (TURB FALSE CACHE BOOL "Enable Turbulence")
set (STATS FALSE CACHE BOOL "Enable Statistics")
set (PLAG FALSE CACHE BOOL "Enable Lagrangian Particles Multiphase Flow")
set (PEUL FALSE CACHE BOOL "Enable Eulerian Multiphase Flow")
set (PERI FALSE CACHE BOOL "Enabled Periodic")
set (INRT FALSE CACHE BOOL "Enable Particle Interactions")
set (RADI FALSE CACHE BOOL "Enable Radiation")
set (SPEC FALSE CACHE BOOL "Enable Species")
add_definitions(-DNO_TECPLOT -DRFLO)
if(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
else()
  find_package(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK -DMPI)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
endif()
#find_library(HDF4_LIB df HINTS /usr/lib/hdf /usr/lib64/hdf)
#find_library(MFHDF_LIB mfhdf HINTS /usr/lib/hdf /usr/lib64/hdf)
find_library(ZLIB z)
find_library(JPEG_LIB jpeg)
find_path(HDF4_INC hdf.h HINTS /usr/include/hdf)

#INCLUDE(CTest)
find_path(ROCCOM_INC com.h HINTS ../../IMPACT/COM/include)
find_path(ROCMAN_INC rocman.h HINTS ../../Rocman/include)

file (GLOB CORE_SRCS ../rocflo/*.F90)
set (LIBFLOFLU_SRCS
  ../libfloflu/AfterUpdateMP.F90
  ../libfloflu/AllocateMemoryWork.F90
  ../libfloflu/AverageVecVar.F90
  ../libfloflu/BcondFarfieldPerf.F90
  ../libfloflu/BcondInflowPerf.F90
  ../libfloflu/BcondInjectionPerf.F90
  ../libfloflu/BcondOutflowPerf.F90
  ../libfloflu/BuildPatchIdentifier.F90
  ../libfloflu/CellCentroid.F90
  ../libfloflu/CellGradientsMP.F90
  ../libfloflu/ControlVolume.F90
  ../libfloflu/ConvectiveFluxes.F90
  ../libfloflu/ConvectiveFluxesMP.F90
  ../libfloflu/DescaleGridSpeeds.F90
  ../libfloflu/ExplicitMultistage.F90
  ../libfloflu/FaceCentroid.F90
  ../libfloflu/FaceVector.F90
  ../libfloflu/GlobalCommunicationMP.F90
  ../libfloflu/InitCommunicationMP.F90
  ../libfloflu/IntegrateSourceTermsMP.F90
  ../libfloflu/MakeNumberedKeys.F90
  ../libfloflu/MixtPerf_C.F90
  ../libfloflu/MixtPerf_D.F90
  ../libfloflu/MixtPerf_E.F90
  ../libfloflu/MixtPerf_G.F90
  ../libfloflu/MixtPerf_H.F90
  ../libfloflu/MixtPerf_M.F90
  ../libfloflu/MixtPerf_P.F90
  ../libfloflu/MixtPerf_R.F90
  ../libfloflu/MixtPerf_T.F90
  ../libfloflu/MixtPerf_Vm.F90
  ../libfloflu/MixtureProperties.F90
  ../libfloflu/NumericalDissipation.F90
  ../libfloflu/NumericalDissipationMP.F90
  ../libfloflu/PerfgasDependentVars.F90
  ../libfloflu/PerfgasGasVars.F90
  ../libfloflu/PerfgasTransportVars.F90
  ../libfloflu/ReadAccelerationSection.F90
  ../libfloflu/ReadFlowSection.F90
  ../libfloflu/ReadForcesSection.F90
  ../libfloflu/ReadFormatsSection.F90
  ../libfloflu/ReadGridMotionSection.F90
  ../libfloflu/ReadInitFlowSection.F90
  ../libfloflu/ReadInputFile.F90
  ../libfloflu/ReadInputUtil.F90
  ../libfloflu/ReadMultigridSection.F90
  ../libfloflu/ReadNumericsSection.F90
  ../libfloflu/ReadPostSection.F90
  ../libfloflu/ReadProbeSection.F90
  ../libfloflu/ReadRandomSection.F90
  ../libfloflu/ReadReferenceSection.F90
  ../libfloflu/ReadThrustSection.F90
  ../libfloflu/ReadTimestepSection.F90
  ../libfloflu/ReadTransformSection.F90
  ../libfloflu/ReadViscositySection.F90
  ../libfloflu/RkInitGeneric.F90
  ../libfloflu/RkInitMP.F90
  ../libfloflu/RkInitPointScalar.F90
  ../libfloflu/RkUpdateGeneric.F90
  ../libfloflu/RkUpdateMP.F90
  ../libfloflu/RkUpdatePointScalar.F90
  ../libfloflu/RungeKuttaMP.F90
  ../libfloflu/ScaleGridSpeeds.F90
  ../libfloflu/ScaleRotateVector.F90
  ../libfloflu/SourceTerms.F90
  ../libfloflu/SourceTermsMP.F90
  ../libfloflu/SplitQuadFace.F90
  ../libfloflu/UpdateBoundaryConditionsMP.F90
  ../libfloflu/UpdateDependentVarsMP.F90
  ../libfloflu/UpdateTbc.F90
  ../libfloflu/UpdateTbcPiecewise.F90
  ../libfloflu/UpdateTbcSinusoidal.F90
  ../libfloflu/UpdateTbcStochastic.F90
  ../libfloflu/UpdateTbcWhitenoise.F90
  ../libfloflu/ViscousFluxes.F90
  ../libfloflu/ViscousFluxesMP.F90
  ../libfloflu/WriteConvergence.F90
  ../libfloflu/WriteProbe.F90
  ../libfloflu/WriteThrust.F90
  ../libfloflu/WriteTotalMass.F90
  ../libfloflu/ZeroDummyCellsMP.F90
  ../libfloflu/ZeroResidualsMP.F90 )
if(STATS)
  set (LIBFLOFLU_SRCS ${LIBFLOFLU_SRCS} ../libfloflu/ReadStatisticSection.F90)
endif()
set (LIBFLO_SRCS
  ../libflo/RFLO_ArcLengthBounds.F90
  ../libflo/RFLO_BoundaryDeformation.F90
  ../libflo/RFLO_CalcCellCentroids.F90
  ../libflo/RFLO_CalcFaceCentroids.F90
  ../libflo/RFLO_CalcGradConnBc.F90
  ../libflo/RFLO_CalcGradDummyConn.F90
  ../libflo/RFLO_CalcGradDummy.F90
  ../libflo/RFLO_CalcGradDummyPhys.F90
  ../libflo/RFLO_CalcGradDummySymm.F90
  ../libflo/RFLO_CalcGradFaces.F90
  ../libflo/RFLO_CalcGradPhysBc.F90
  ../libflo/RFLO_CalcGradVector.F90
  ../libflo/RFLO_ChangeInteriorGrid.F90
  ../libflo/RFLO_CheckValidity.F90
  ../libflo/RFLO_CopyBoundaryData.F90
  ../libflo/RFLO_CopyEdgeFaceNorm.F90
  ../libflo/RFLO_CopyEdgeFaceParal.F90
  ../libflo/RFLO_CopyGeometryDummy.F90
  ../libflo/RFLO_CopyTopologyLevels.F90
  ../libflo/RFLO_DerivedInputValues.F90
  ../libflo/RFLO_EdgeDeformation.F90
  ../libflo/RFLO_EdgeDeformationStraight.F90
  ../libflo/RFLO_ExchangeDnodeCopy.F90
  ../libflo/RFLO_ExchangeDnodeSend.F90
  ../libflo/RFLO_ExchangeDnodeRecv.F90
  ../libflo/RFLO_ExtrapIntCellVec.F90
  ../libflo/RFLO_ExtrapIntCellScal.F90
  ../libflo/RFLO_GenerateCoarseGrids.F90
  ../libflo/RFLO_GetDimensions.F90
  ../libflo/RFLO_GetPatchDimens.F90
  ../libflo/RFLO_GetPatchMapping.F90
  ../libflo/RFLO_GridRemesh.F90
  ../libflo/RFLO_InitInputValues.F90
  ../libflo/RFLO_InterpolDistrib.F90
  ../libflo/RFLO_RandomInit.F90
  ../libflo/RFLO_ReadBcFarfSection.F90
  ../libflo/RFLO_ReadBcFromFile.F90
  ../libflo/RFLO_ReadBcInflowTotAngSection.F90
  ../libflo/RFLO_ReadBcInflowVelSection.F90
  ../libflo/RFLO_ReadBcInjectMrateSection.F90
  ../libflo/RFLO_ReadBcInjectAPNSection.F90
  ../libflo/RFLO_ReadBcInputFile.F90
  ../libflo/RFLO_ReadBcNoslipSection.F90
  ../libflo/RFLO_ReadBcOutflowSection.F90
  ../libflo/RFLO_ReadBcSlipWallSection.F90
  ../libflo/RFLO_ReadDataFile.F90
  ../libflo/RFLO_ReadGrid.F90
  ../libflo/RFLO_ReadGridRegion.F90
  ../libflo/RFLO_ReadRandomState.F90
  ../libflo/RFLO_ReadRegionTopology.F90
  ../libflo/RFLO_ReadSolution.F90
  ../libflo/RFLO_ReadSolutionRegion.F90
  ../libflo/RFLO_ReadTbcInputFile.F90
  ../libflo/RFLO_ReadTbcSection.F90
  ../libflo/RFLO_SetMstageCoeffs.F90
  ../libflo/RFLO_Tfint.F90
  ../libflo/RFLO_WriteDataFile.F90
  ../libflo/RFLO_WriteGrid.F90
  ../libflo/RFLO_WriteGridRegion.F90
  ../libflo/RFLO_WriteRandomState.F90
  ../libflo/RFLO_WriteRegionTopology.F90
  ../libflo/RFLO_WriteSolution.F90
  ../libflo/RFLO_WriteSolutionRegion.F90
  ../libflo/RFLO_ZeroDummyCells.F90)
if (STATS)
  set (LIBFLO_SRCS ${LIBFLO_SRCS} ../libflo/RFLO_ReadStat.F90
    ../libflo/RFLO_ReadStatRegion.F90
    ../libflo/RFLO_WriteStat.F90)
endif()

set ( MODFLOFLU_SRCS
  ../modfloflu/ModDataTypes.F90
  ../modfloflu/ModError.F90
  ../modfloflu/ModGlobal.F90
  ../modfloflu/ModMPI.F90
  ../modfloflu/ModParameters.F90
  ../modfloflu/ModDataStruct.F90
  ../modfloflu/ModBndPatch.F90
  ../modfloflu/ModBuildFileNames.F90
  ../modfloflu/ModGrid.F90
  ../modfloflu/ModInteract.F90
  ../modfloflu/ModInterfacesBcond.F90
  ../modfloflu/ModInterfacesIO.F90
  ../modfloflu/ModInterfacesMixt.F90
  ../modfloflu/ModInterfacesUtil.F90
  ../modfloflu/ModInterfaces.F90
  ../modfloflu/ModMaterials.F90
  ../modfloflu/ModMixture.F90
  ../modfloflu/ModPartEul.F90
  ../modfloflu/ModPartLag.F90
  ../modfloflu/ModPeriodic.F90
  ../modfloflu/ModRadiation.F90
  ../modfloflu/ModRandom.F90
  ../modfloflu/ModSortSearch.F90
  ../modfloflu/ModSpecies.F90
  ../modfloflu/ModStencil.F90
  ../modfloflu/ModTools.F90
  ../modfloflu/ModTurbulence.F90 )

if(ROCSTAR)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModRocstar.F90 )
endif()
if(INRT)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesInteract.F90 )
endif()
if(PEUL)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesEulerian.F90 )
endif()
if(PLAG)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesLagrangian.F90 )
endif()
if(PERI)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesPeriodic.F90 )
endif()
if(RADI)
  list (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesRadiation.F90 )
endif()
IF(SPEC)
  LIST (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesSpecies.F90 )
ENDIF()
IF(STATS)
  LIST (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesStatistics.F90
    ../modfloflu/ModStatsRoutines.F90 )
ENDIF()
IF(TURB)
  LIST (APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesTurbulence.F90 )
ENDIF()

SET (MODFLO_SRCS
  ../modflo/ModIndexing.F90
  ../modflo/RFLO_ModInterfacesExternal.F90
  ../modflo/RFLO_ModInterfacesLibrary.F90
  ../modflo/RFLO_ModInterfacesSolver.F90
  ../modflo/RFLO_ModInitialization.F90
  ../modflo/RFLO_ModBoundaryConditions.F90
  ../modflo/RFLO_ModDegenerateCornEdge.F90
  ../modflo/RFLO_ModElliptSmoothing.F90
  ../modflo/RFLO_ModExtrapolation.F90
  ../modflo/RFLO_ModFiniteDifference.F90
  ../modflo/RFLO_ModForcesMoments.F90
  ../modflo/RFLO_ModGridControlMap.F90
  ../modflo/RFLO_ModGridMetrics.F90
  ../modflo/RFLO_ModGridRegionShape.F90
  ../modflo/RFLO_ModLaplaceSmoothing.F90
  ../modflo/RFLO_ModMoveGridFrame.F90
  ../modflo/RFLO_ModMoveGridElliptGlo.F90
  ../modflo/RFLO_ModMoveGridElliptFra.F90
  ../modflo/RFLO_ModMoveGridUtil.F90
  ../modflo/RFLO_ModPatchAeroCoeffs.F90
  ../modflo/RFLO_ModRestartInfo.F90
  ../modflo/RFLO_ModVectorTensor.F90
  ../modflo/RFLO_ModVolMeshSmoothing.F90 )
IF(STATS)
  LIST (APPEND MODFLO_SRCS ../modflo/RFLO_ModStatsBoundaryConditions.F90 )
ENDIF()


SET (TURB_SRCS
  ../rocturb/TURB_ModInterfaces.F90
  ../rocturb/TURB_ModParameters.F90
  ../rocturb/TURB_AllocateMemory.F90
  ../rocturb/TURB_BuildVersionString.F90
  ../rocturb/TURB_CalcMetrics.F90
  ../rocturb/TURB_CalcStrainRate.F90
  ../rocturb/TURB_CalcVortic.F90
  ../rocturb/TURB_CheckParamInput.F90
  ../rocturb/TURB_DerivedInputValues.F90
  ../rocturb/TURB_EmsInit.F90
  ../rocturb/TURB_GetModelStressCell.F90
  ../rocturb/TURB_GetTvCell.F90
  ../rocturb/TURB_InitInputValues.F90
  ../rocturb/TURB_InitSolution.F90
  ../rocturb/TURB_LesCalcEddyVis.F90
  ../rocturb/TURB_LesCoefDynMixd.F90
  ../rocturb/TURB_LesCoefDynSmag.F90
  ../rocturb/TURB_LesCommunication.F90
  ../rocturb/TURB_LesContract.F90
  ../rocturb/TURB_LesEsgModel4.F90
  ../rocturb/TURB_LesFluxFixSmag.F90
  ../rocturb/TURB_LesFluxScalSim.F90
  ../rocturb/TURB_LesGetEddyVis.F90
  ../rocturb/TURB_LesHij.F90
  ../rocturb/TURB_LesLij.F90
  ../rocturb/TURB_LesMij.F90
  ../rocturb/TURB_LesRkInit.F90
  ../rocturb/TURB_LesTestRhoV.F90
  ../rocturb/TURB_PrintUserInput.F90
  ../rocturb/TURB_RansConvectiveFluxes.F90
  ../rocturb/TURB_RansEmsInit.F90
  ../rocturb/TURB_RansEmsUpdate.F90
  ../rocturb/TURB_RansNumericalDissipation.F90
  ../rocturb/TURB_RansRkInit.F90
  ../rocturb/TURB_RansSAGetEddyVis.F90
  ../rocturb/TURB_RansSASourceTerms.F90
  ../rocturb/TURB_RansSAVisFlux.F90
  ../rocturb/TURB_RansSAVisFluxPatch.F90
  ../rocturb/TURB_RansSourceTerms.F90
  ../rocturb/TURB_RansTotalTv.F90
  ../rocturb/TURB_RansWallDistOVPatch.F90
  ../rocturb/TURB_RansZeroDummyCells.F90
  ../rocturb/TURB_ReadBcInputFile.F90
  ../rocturb/TURB_ReadInputFile.F90
  ../rocturb/TURB_ReadTurbSection.F90
  ../rocturb/TURB_RkInit.F90
  ../rocturb/TURB_SolutionUpdate.F90
  ../rocturb/TURB_UserInput.F90
  ../rocturb/TURB_VFluxHybrid.F90
  ../rocturb/TURB_VFluxHybridPatch.F90
  ../rocturb/TURB_VisFluxEddy.F90
  ../rocturb/TURB_VisFluxEddyPatch.F90
  ../rocturb/TURB_WlmFluxPatch.F90
  ../rocturb/TURB_WlmInitia.F90
  ../rocturb/TURB_WlmReyAnalogy.F90
  ../rocturb/TURB_WlmTauWallMapping.F90
  ../rocturb/TURB_WlmUpdate.F90
  ../rocturb/TURB_WlmUpdateBndlay.F90
  ../rocturb/TURB_coRansWallDistOVFlo.F90
  ../rocturb/TURB_coViscousFluxesFlo.F90
  ../rocturb/TURB_coWlmReadBcSectionFlo.F90
  ../rocturb/TURB_floExtrapIntCellScal.F90
  ../rocturb/TURB_floExtrapIntCellVec.F90
  ../rocturb/TURB_floExtrapolCellVec.F90
  ../rocturb/TURB_floExtrapolFaceVec.F90
  ../rocturb/TURB_floFaceVolume.F90
  ../rocturb/TURB_floFaceWidthDummyConn.F90
  ../rocturb/TURB_floFaceWidthDummy.F90
  ../rocturb/TURB_floFaceWidthDummyPhys.F90
  ../rocturb/TURB_floFaceWidth.F90
  ../rocturb/TURB_floLesAverageFace.F90
  ../rocturb/TURB_floLesGenC2F.F90
  ../rocturb/TURB_floLesGenCoCC.F90
  ../rocturb/TURB_floLesGenCoCCUtil.F90
  ../rocturb/TURB_floLesGenCoFCUtil.F90
  ../rocturb/TURB_floLesGenCoFF.F90
  ../rocturb/TURB_floLesGenCoFFUtil.F90
  ../rocturb/TURB_floLesGenFiltCC.F90
  ../rocturb/TURB_floLesGenFiltCCI.F90
  ../rocturb/TURB_floLesGenFiltCCJ.F90
  ../rocturb/TURB_floLesGenFiltCCK.F90
  ../rocturb/TURB_floLesGenFiltFF.F90
  ../rocturb/TURB_floLesGenFiltFFI.F90
  ../rocturb/TURB_floLesGenFiltFFJ.F90
  ../rocturb/TURB_floLesGenFiltFFK.F90
  ../rocturb/TURB_floLesUniFiltCC.F90
  ../rocturb/TURB_floLesUniFiltCCI.F90
  ../rocturb/TURB_floLesUniFiltCCJ.F90
  ../rocturb/TURB_floLesUniFiltCCK.F90
  ../rocturb/TURB_floLesUniFiltFF.F90
  ../rocturb/TURB_floLesUniFiltFFI.F90
  ../rocturb/TURB_floLesUniFiltFFJ.F90
  ../rocturb/TURB_floLesUniFiltFFK.F90
  ../rocturb/TURB_floRansBcondInflow.F90
  ../rocturb/TURB_floRansBcondInjection.F90
  ../rocturb/TURB_floRansBcondNoslipWall.F90
  ../rocturb/TURB_floRansBcondSymmetry.F90
  ../rocturb/TURB_floRansBcondZeroGrad.F90
  ../rocturb/TURB_floRansCentralDissipation.F90
  ../rocturb/TURB_floRansCorrCornEdgeCells.F90
  ../rocturb/TURB_floRansExchangeDummyConf.F90
  ../rocturb/TURB_floRansExchCornEdgeCells.F90
  ../rocturb/TURB_floRansRecvDummyVals.F90
  ../rocturb/TURB_floRansRecvCornEdgeCells.F90
  ../rocturb/TURB_floRansSACentFlux.F90
  ../rocturb/TURB_floRansSACentFluxPatch.F90
  ../rocturb/TURB_floRansSARoe1stFlux.F90
  ../rocturb/TURB_floRansSARoe2ndFlux.F90
  ../rocturb/TURB_floRansSARoeFluxPatch.F90
  ../rocturb/TURB_floRansSendDummyConf.F90
  ../rocturb/TURB_floRansSendCornEdgeCells.F90
  ../rocturb/TURB_floRansSetCornEdgeCells.F90
  ../rocturb/TURB_floWlmMetric.F90
  ../rocturb/TURB_floWlmUpdateLoglay.F90
  ../rocturb/TURB_rFLO_RansAllocDataBuffers.F90
  ../rocturb/TURB_rFLO_RansBndConditionsRecv.F90
  ../rocturb/TURB_rFLO_RansBndConditionsSend.F90
  ../rocturb/TURB_rFLO_RansBndConditionsSet.F90
  ../rocturb/TURB_rFLO_RansClearSendRequests.F90
  ../rocturb/TURB_rFLO_RansResSmoothing.F90
  ../rocturb/TURB_rFLO_RansResSmoothingCoeff.F90
  ../rocturb/TURB_rFLO_RansSpectralRadii.F90
  ../rocturb/TURB_rFLO_ReadSolution.F90
  ../rocturb/TURB_rFLO_ReadSolutionRegion.F90
  ../rocturb/TURB_rFLO_WriteSolution.F90 )

IF(STATS)
  LIST (APPEND TURB_SRCS   ../rocturb/TURB_StatCCollector.F90
    ../rocturb/TURB_StatFCollector.F90
    ../rocturb/TURB_StatMapping.F90 )
ENDIF()
FILE(GLOB RADI_SRCS ../rocrad/*.F90)
SET(PLAG_SRCS
  ../rocpart/PLAG_BuildVersionString.F90
  ../rocpart/PLAG_ModInterfaces.F90
  ../rocpart/PLAG_ModInterfacesCECells.F90
  ../rocpart/PLAG_ModParameters.F90
  ../rocpart/PLAG_ModCheckVars.F90
  ../rocpart/PLAG_ModRkInit.F90
  ../rocpart/PLAG_AllocateDataBuffers.F90
  ../rocpart/PLAG_AllocateMemory.F90
  ../rocpart/PLAG_AllocateMemoryTile.F90
  ../rocpart/PLAG_AppendDataFromBuffers.F90
  ../rocpart/PLAG_BoundaryConditionsSet.F90
  ../rocpart/PLAG_BufferDataRecv.F90
  ../rocpart/PLAG_BufferDataSend.F90
  ../rocpart/PLAG_BufferSizeRecv.F90
  ../rocpart/PLAG_BufferSizeSend.F90
  ../rocpart/PLAG_CalcBreakup.F90
  ../rocpart/PLAG_CalcDerivedVariables.F90
  ../rocpart/PLAG_CalcFaceCentroids.F90
  ../rocpart/PLAG_CalcRhsPosition.F90
  ../rocpart/PLAG_CECellsAllocateData.F90
  ../rocpart/PLAG_CECellsClearRequestsData.F90
  ../rocpart/PLAG_CECellsClearRequestsSize.F90
  ../rocpart/PLAG_CECellsDeallocateData.F90
  ../rocpart/PLAG_CECellsExchange.F90
  ../rocpart/PLAG_CECellsFaceCentroids.F90
  ../rocpart/PLAG_CECellsFaceVectors.F90
  ../rocpart/PLAG_CECellsGetBufferSize.F90
  ../rocpart/PLAG_CECellsLoadDataWrapper.F90
  ../rocpart/PLAG_CECellsRecvData.F90
  ../rocpart/PLAG_CECellsRecvSize.F90
  ../rocpart/PLAG_CECellsSendData.F90
  ../rocpart/PLAG_CECellsSendRecvWrapper.F90
  ../rocpart/PLAG_CECellsSendSize.F90
  ../rocpart/PLAG_CECellsWrapper.F90
  ../rocpart/PLAG_CheckUserInput.F90
  ../rocpart/PLAG_ClearDataSendRequests.F90
  ../rocpart/PLAG_ClearSizeSendRequests.F90
  ../rocpart/PLAG_CopyFaceVectors.F90
  ../rocpart/PLAG_CornCellsLoadData.F90
  ../rocpart/PLAG_CornCellsLoadSendBuff.F90
  ../rocpart/PLAG_DeallocateDataBuffers.F90
  ../rocpart/PLAG_DerivedInputValues.F90
  ../rocpart/PLAG_EdgeCellsLoadData.F90
  ../rocpart/PLAG_EdgeCellsLoadSendBuff.F90
  ../rocpart/PLAG_GetCellIndices.F90
  ../rocpart/PLAG_GetCellIndicesOutflow.F90
  ../rocpart/PLAG_InCellTest.F90
  ../rocpart/PLAG_InCellTestRobust.F90
  ../rocpart/PLAG_InitInputValues.F90
  ../rocpart/PLAG_InitMemory.F90
  ../rocpart/PLAG_InitSolution.F90
  ../rocpart/PLAG_InjcEjectParticle.F90
  ../rocpart/PLAG_InjcMakeParticle.F90
  ../rocpart/PLAG_InjcSetInjection.F90
  ../rocpart/PLAG_InjcTileCalcRhs.F90
  ../rocpart/PLAG_InjcTileInitialize.F90
  ../rocpart/PLAG_InjcTileRKUpdate.F90
  ../rocpart/PLAG_InjcTileUpdate.F90
  ../rocpart/PLAG_InjcTileZeroRhs.F90
  ../rocpart/PLAG_IntrpMixtProperties.F90
  ../rocpart/PLAG_NonCvUpdate.F90
  ../rocpart/PLAG_PatchBufferSendRecv.F90
  ../rocpart/PLAG_PatchExchangeConf.F90
  ../rocpart/PLAG_PatchGetBufferSize.F90
  ../rocpart/PLAG_PatchLoadDataBuffers.F90
  ../rocpart/PLAG_PatchRemoveDataOutflow.F90
  ../rocpart/PLAG_PatchUpdate.F90
  ../rocpart/PLAG_PatchUpdateWrapper.F90
  ../rocpart/PLAG_PrintUserInput.F90
  ../rocpart/PLAG_ReadDisPartInitSection.F90
  ../rocpart/PLAG_ReadDisPartnContSection.F90
  ../rocpart/PLAG_ReadDisPartSection.F90
  ../rocpart/PLAG_ReadInputFile.F90
  ../rocpart/PLAG_ReadPdfFromFile.F90
  ../rocpart/PLAG_ReadSolution.F90
  ../rocpart/PLAG_RkInit.F90
  ../rocpart/PLAG_RkUpdateWrapper.F90
  ../rocpart/PLAG_Update.F90
  ../rocpart/PLAG_UpdateDataStruct.F90
  ../rocpart/PLAG_UserInput.F90
  ../rocpart/PLAG_WallBounce.F90
  ../rocpart/PLAG_WriteSolution.F90
  ../rocpart/PLAG_ZeroRhs.F90
  ../rocpart/PLAG_RFLO_SetMetrics.F90
  ../rocpart/PLAG_RFLO_ClearSendRequests.F90
  ../rocpart/PLAG_RFLO_FindGridMapping.F90
  ../rocpart/PLAG_RFLO_RecvMetrics.F90
  ../rocpart/PLAG_RFLO_SendMetrics.F90
  ../rocpart/PLAG_RFLO_SendRecvMetrics.F90)

IF(STATS)
  LIST(APPEND PLAG_SRCS ../rocpart/PLAG_ModEulerian.F90
    ../rocpart/PLAG_ModStats.F90
    ../rocpart/PLAG_RFLO_ModStats.F90
    ../rocpart/PLAG_StatMapping.F90)
ENDIF()
FILE(GLOB PEUL_SRCS ../rocsmoke/*.F90)
FILE(GLOB INRT_SRCS ../rocinteract/*.F90)
SET(ROCSTAR_SRCS
  ../rocstar/RFLO_FlowSolverDummy.F90
  ../rocstar/RFLO_GetBoundaryValues.F90
  ../rocstar/RFLO_GetDeformation.F90
  ../rocstar/RFLO_InitRocstarInterface.F90
  ../rocstar/RFLO_SendBoundaryValues.F90
  ../rocstar/RFLO_SendBoundaryValuesAlpha.F90
  ../rocstar/RFLO_UpdateInbuffGm.F90
  ../rocstar/Fluid_finalize.F90
  ../rocstar/Fluid_preHdfOutput.F90
  ../rocstar/Fluid_postHdfOutput.F90
  ../rocstar/Rocflo_load_module.F90
  ../rocstar/Rocflo_unload_module.F90)
IF(PEUL)
  LIST(APPEND ROCSTAR_SRCS ../rocstar/PEUL_InitRocstarInterface.F90)
ENDIF()
if(PLAG)
  LIST(APPEND ROCSTAR_SRCS ../rocstarPLAG_InitRocstarInterface.F90
    ../rocstar/PLAG_SetSizeRocstar.F90)
endif()
if(RADI)
  LIST(APPEND ROCSTAR_SRCS ../rocstar/RADI_InitRocstarInterface.F90 )
endif()
if(TURB)
  LIST(APPEND ROCSTAR_SRCS ../rocstar/TURB_InitRocstarInterface.F90)
endif()
LIST(APPEND ROCSTAR_SRCS ../rocstar/RandInitRocstarInterface.F90)
if(STATS)
  LIST(APPEND ROCSTAR_SRCS ../rocstar/STAT_RocstarStatNaming.F90)
endif()
IF(STATS)
  LIST(APPEND ROCSTAR_SRCS ../rocstar/STAT_RocstarStatNaming.F90)
ENDIF()
SET(FLOPREP_SRCS
  ../utilities/rocflo/prep/PREP_ModInterfaces.F90
  ../utilities/rocflo/prep/PREP_ModParameters.F90
  ../utilities/rocflo/prep/PREP_ModBcDistribution.F90
  ../utilities/rocflo/prep/PREP_BuildVersionString.F90
  ../utilities/rocflo/prep/PREP_CheckBcValidity.F90
  ../utilities/rocflo/prep/PREP_GetGrid.F90
  ../utilities/rocflo/prep/PREP_InitializeFlowField.F90
  ../utilities/rocflo/prep/PREP_InitInputValues.F90
  ../utilities/rocflo/prep/PREP_Main.F90
  ../utilities/rocflo/prep/PREP_PrintPrepInput.F90
  ../utilities/rocflo/prep/PREP_ReadInitflowSection.F90
  ../utilities/rocflo/prep/PREP_ReadInputFile.F90
  ../utilities/rocflo/prep/PREP_ReadBcInputFile.F90)
IF(ROCSTAR)
  LIST(APPEND FLOPREP_SRCS ../utilities/rocflo/prep/PREP_RocstarInitSolution.F90
    ../utilities/rocflo/prep/PREP_RocstarWriteRocinout.F90
    ../utilities/rocflo/prep/PREP_RocstarWriteSolution.F90
    ../utilities/rocflo/prep/PREP_InitRocstarInterface.F90)
ENDIF()

SET(LIB_SRCS ${MODFLO_SRCS} ${MODFLOFLU_SRCS} ${LIBFLO_SRCS} ${LIBFLOFLU_SRCS}
  ${CORE_SRCS})
IF(ROCSTAR)
  add_definitions(-DROCSTAR -DGENX)
  LIST(APPEND LIB_SRCS ${ROCSTAR_SRCS})
ENDIF()
IF(RADI)
  add_definitions(-DRADI)
  LIST(APPEND LIB_SRCS ${RADI_SRCS})
ENDIF()
IF(TURB)
  add_definitions(-DTURB)
  LIST(APPEND LIB_SRCS ${TURB_SRCS})
ENDIF()
IF(PLAG)
  add_definitions(-DPLAG)
  LIST(APPEND LIB_SRCS ${PLAG_SRCS})
ENDIF()
if(PEUL)
  add_definitions(-DPEUL)
  list(APPEND LIB_SRCS ${PEUL_SRCS})
endif()
if(INRT)
  add_definitions(-DINRT)
  list(APPEND LIB_SRCS ${INRT_SRCS})
endif()
if(STATS)
  add_definitions(-DSTATS)
endif()
if(PERI)
  add_definitions(-DPERI)
endif()
if(SPEC)
  add_definitions(-DSPEC)
endif()
set_source_files_properties(${LIB_SRCS} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
#
#UTIL_OBJ=hash_table.o read.a write.a util.o
#OBJ= adj.o face.o buildface.o split.o volume.o
#	readinp.o writeflo.o writetop.o writemblock.o
#	$(UTIL_OBJ)
#HEADERS=adj.h face.h util.h hash_table.h vector3d.h gridutil.h
#DEST=$(BINPRE)makeflo    # Prepend $(BINPRE) to put into bin directory
#RTST= $(UTIL_OBJ)
#ATST=read.o $(UTIL_OBJ)
#
#CPPFLAGS+=-DUSE_HDF=1 $(HDFINC)
#LFLAGS=$(LDFLAGS) $(HDFLIB) $(SYSLIBS) $(DEPLIBS)
#

set(MAKEFLO_DIR ../utilities/rocflo/makeflo)
file(GLOB MAKEFLO_CPP_FILES ${MAKEFLO_DIR}/*.cpp)
file(GLOB MAKEFLO_F90_FILES ${MAKEFLO_DIR}/*.f90)
set_source_files_properties(${MAKEFLO_CPP_FILES} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_source_files_properties(${MAKEFLO_F90_FILES} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set(MAKEFLO_READ_SRCS ${MAKEFLO_DIR}/read.cpp ${MAKEFLO_DIR}/readgrd.cpp ${MAKEFLO_DIR}/readreg.cpp
                      ${MAKEFLO_DIR}/readmsh.f90 ${MAKEFLO_DIR}/readhdf.cpp)
set(MAKEFLO_WRITE_SRCS ${MAKEFLO_DIR}/write.cpp ${MAKEFLO_DIR}/writemsh.f90 ${MAKEFLO_DIR}/writehdf.cpp
                       ${MAKEFLO_DIR}/writegrda.f90 ${MAKEFLO_DIR}/writegrdb.f90)
set(MAKEFLO_UTIL_SRCS ${MAKEFLO_DIR}/hash_table.cpp ${MAKEFLO_DIR}/util.cpp)
set(MAKEFLO_CORE_SRCS ${MAKEFLO_DIR}/adj.cpp ${MAKEFLO_DIR}/face.cpp ${MAKEFLO_DIR}/buildface.cpp
                      ${MAKEFLO_DIR}/split.cpp ${MAKEFLO_DIR}/volume.cpp ${MAKEFLO_DIR}/readinp.cpp
                      ${MAKEFLO_DIR}/writeflo.cpp ${MAKEFLO_DIR}/writetop.cpp ${MAKEFLO_DIR}/writemblock.cpp
                      ${MAKEFLO_UTIL_SRCS})
set(MAKEFLO_SRCS ${MAKEFLO_CPP_FILES} ${MAKEFLO_F90_FILES})
set(UTIL_SRCS ${MAKEFLO_SRCS})
set(ALL_SRCS ${LIB_SRCS} ${UTIL_SRCS})
set_source_files_properties(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC" )
#set(TEST_SRCS src/MANTest.C)

# rpath settings
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories( ../include ${ROCCOM_INC} ${ROCMAN_INC} ${HDF4_INC})
add_definitions(-DUSE_HDF=1)

if(NOT BUILD_STATIC)
  add_library(Rocflo SHARED ${LIB_SRCS})
else()
  add_library(Rocflo STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
endif()
target_link_libraries(Rocflo ${IMPACT_LIB} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

#add_library(makeflofio ${MAKEFLO_DIR}/readmsh.f90 ${MAKEFLO_DIR}/writemsh.f90 ${MAKEFLO_DIR}/writegrda.f90 ${MAKEFLO_DIR}/writegrdb.f90)
#add_library(read ${MAKEFLO_DIR}/read.cpp ${MAKEFLO_DIR}/readgrd.cpp ${MAKEFLO_DIR}/readreg.cpp ${MAKEFLO_DIR}/readhdf.cpp)
#target_link_libraries(read ${CMAKE_Fortran_STANDARD_LIBRARIES})
#add_library(write ${MAKEFLO_DIR}/write.cpp ${MAKEFLO_DIR}/writehdf.cpp)
#target_link_libraries(write ${CMAKE_Fortran_STANDARD_LIBRARIES})
add_executable(makeflo ${MAKEFLO_DIR}/makeflo.cpp ${MAKEFLO_CORE_SRCS} ${MAKEFLO_DIR}/read.cpp ${MAKEFLO_DIR}/readgrd.cpp ${MAKEFLO_DIR}/readreg.cpp ${MAKEFLO_DIR}/readhdf.cpp ${MAKEFLO_DIR}/write.cpp ${MAKEFLO_DIR}/writehdf.cpp ${MAKEFLO_DIR}/readmsh.f90 ${MAKEFLO_DIR}/writemsh.f90 ${MAKEFLO_DIR}/writegrda.f90 ${MAKEFLO_DIR}/writegrdb.f90)
target_link_libraries(makeflo ${MFHDF_LIB} ${HDF4_LIB} ${JPEG_LIB} ${ZLIB} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
add_executable(makemblock ${MAKEFLO_DIR}/makemblock.cpp ${MAKEFLO_CORE_SRCS} ${MAKEFLO_DIR}/read.cpp ${MAKEFLO_DIR}/readgrd.cpp ${MAKEFLO_DIR}/readreg.cpp ${MAKEFLO_DIR}/readhdf.cpp ${MAKEFLO_DIR}/write.cpp ${MAKEFLO_DIR}/writehdf.cpp ${MAKEFLO_DIR}/readmsh.f90 ${MAKEFLO_DIR}/writemsh.f90 ${MAKEFLO_DIR}/writegrda.f90 ${MAKEFLO_DIR}/writegrdb.f90)
target_link_libraries(makemblock ${MFHDF_LIB} ${HDF4_LIB} ${JPEG_LIB} ${ZLIB} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
add_executable(simpleblock ${MAKEFLO_DIR}/simpleblock.cpp ${MAKEFLO_CORE_SRCS} ${MAKEFLO_DIR}/read.cpp ${MAKEFLO_DIR}/readgrd.cpp ${MAKEFLO_DIR}/readreg.cpp ${MAKEFLO_DIR}/readhdf.cpp ${MAKEFLO_DIR}/write.cpp ${MAKEFLO_DIR}/writehdf.cpp ${MAKEFLO_DIR}/readmsh.f90 ${MAKEFLO_DIR}/writemsh.f90 ${MAKEFLO_DIR}/writegrda.f90 ${MAKEFLO_DIR}/writegrdb.f90)
target_link_libraries(simpleblock ${MFHDF_LIB} ${HDF4_LIB} ${JPEG_LIB} ${ZLIB} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
add_executable(rfloprep ${FLOPREP_SRCS})
set_property(TARGET rfloprep PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(rfloprep Rocflo ${IMPACT_LIB} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

set_target_properties(Rocflo PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
set_target_properties(makemblock PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
set_target_properties(simpleblock PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
set_target_properties(rfloprep PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
set_target_properties(makeflo PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")

# Unit test programs
#read.tst: $(RTST) read.cpp
#	$(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c read.cpp
#	$(LD) $(OPTS) read.o $(RTST) -o $@ $(LFLAGS)

#adj.tst: $(ATST) adj.cpp
#	$(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c adj.cpp
#	$(LD) $(OPTS) adj.o $(ATST) -o $@ $(LFLAGS)

#ADD_TEST(RunAllTests ${EXECUTABLE_OUTPUT_PATH}/testx testresults.txt)
#ADD_TEST(TestObject:ValidHandle ${TEST_RESULTS} TestObject:ValidHandle testresults.txt)
#ADD_TEST(TestObject:Exists ${TEST_RESULTS} TestObject:Exists testresults.txt)

install(TARGETS makemblock simpleblock rfloprep Rocflo makeflo RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
