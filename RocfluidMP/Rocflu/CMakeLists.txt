cmake_minimum_required(VERSION 3.0)
# project(Rocflu C CXX Fortran)

option(RFLU_ROCSTAR "Build as Rocstar Module" ON)
option(RFLU_TURB "Enable Turbulence" ON)
option(RFLU_STATS "Enable Statistics" ON)
option(RFLU_INRT "Enable Particle Interactions" ON)
#[[
option(RFLU_PLAG "Enable Lagrangian Particles Multiphase Flow" OFF)  # Incomplete
option(RFLU_PEUL "Enable Eulerian Multiphase Flow" OFF)  # Missing source
option(RFLU_PERI "Enable Periodic" OFF)  # Missing source
option(RFLU_RADI "Enable Radiation" OFF)  # Missing source
option(RFLU_SPEC "Enable Species" OFF)  # Missing source
option(PETSC "Enable PETSc" OFF)
]]

add_definitions(-DNO_TECPLOT -DRFLU)

# MPI Found at the top level CMake
if(ENABLE_MPI)
#  find_package(MPI REQUIRED)
  add_definitions(-DMPI)
endif(ENABLE_MPI)

find_package(LAPACK REQUIRED)
find_library(METIS_LIB metis)
##include(CTest)

file(GLOB CORE_SRCS ../rocflu/*.F90)

set(LIBFLOFLU_SRCS
    ../libfloflu/AfterUpdateMP.F90
    ../libfloflu/BcondFarfieldPerf.F90
    ../libfloflu/BcondInflowPerf.F90
    ../libfloflu/BcondInflowPerf_GL.F90
    ../libfloflu/BcondInjectionPerf.F90
    ../libfloflu/BcondOutflowPerf.F90
    ../libfloflu/BcondOutflowPerf_GL.F90
    ../libfloflu/BuildPatchIdentifier.F90
    ../libfloflu/CellCentroid.F90
    ../libfloflu/CellGradientsMP.F90
    ../libfloflu/ControlVolume.F90
    ../libfloflu/ConvectiveFluxes.F90
    ../libfloflu/ConvectiveFluxesMP.F90
    ../libfloflu/DescaleGridSpeeds.F90
    ../libfloflu/FaceCentroid.F90
    ../libfloflu/FaceVector.F90
    ../libfloflu/IntegrateSourceTermsMP.F90
    ../libfloflu/InitCommunicationMP.F90
    ../libfloflu/MakeNumberedKeys.F90
    ../libfloflu/MixtGasLiq_C.F90
    ../libfloflu/MixtGasLiq_E.F90
    ../libfloflu/MixtGasLiq_P.F90
    ../libfloflu/MixtLiq_C.F90
    ../libfloflu/MixtLiq_D.F90
    ../libfloflu/MixtPerf_C.F90
    ../libfloflu/MixtPerf_Cv.F90
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
    ../libfloflu/ReadMixtureSection.F90
    ../libfloflu/ReadMultigridSection.F90
    ../libfloflu/ReadNumericsSection.F90
    ../libfloflu/ReadPostSection.F90
    ../libfloflu/ReadPrepSection.F90
    ../libfloflu/ReadProbeSection.F90
    ../libfloflu/ReadRandomSection.F90
    ../libfloflu/ReadReferenceSection.F90
    ../libfloflu/ReadRocketSection.F90
    ../libfloflu/ReadThrustSection.F90
    ../libfloflu/ReadTimestepSection.F90
    ../libfloflu/ReadTimeZoomingSection.F90
    ../libfloflu/ReadTransformSection.F90
    ../libfloflu/ReadViscositySection.F90
    ../libfloflu/ReflectPosition.F90
    ../libfloflu/ReflectVector.F90
    ../libfloflu/RkInitGeneric.F90
    ../libfloflu/RkInitMP.F90
    ../libfloflu/RkInitSD.F90
    ../libfloflu/RkUpdateGeneric.F90
    ../libfloflu/RkUpdateMP.F90
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
    ../libfloflu/WriteTotalMass.F90
    ../libfloflu/ZeroDummyCellsMP.F90
    ../libfloflu/ZeroResidualsMP.F90)
if(RFLU_STATS)
  list(APPEND LIBFLOFLU_SRCS ../libfloflu/ReadStatisticSection.F90)
endif()

set(LIBFLU_SRCS
    ../libflu/clop.C
    ../libflu/RFLU_AllocateMemoryTbc.F90
    ../libflu/RFLU_BuildDataStruct.F90
    ../libflu/RFLU_CheckDerivedUserInput.F90
    ../libflu/RFLU_CheckPositivity.F90
    ../libflu/RFLU_CheckPositivity_GL.F90
    ../libflu/RFLU_CheckPositivityWrapper.F90
    ../libflu/RFLU_CheckUserInput.F90
    ../libflu/RFLU_CheckValidity.F90
    ../libflu/RFLU_CheckValidity_GL.F90
    ../libflu/RFLU_CheckValidityWrapper.F90
    ../libflu/RFLU_ClosePostInfo.F90
    ../libflu/RFLU_CloseRestartInfo.F90
    ../libflu/RFLU_ComputeDCUHREInfo.F90
    ../libflu/RFLU_CreateGrid.F90
    ../libflu/RFLU_CreateRegions.F90
    ../libflu/RFLU_DecideNeedBGradFace.F90
    ../libflu/RFLU_DecidePrint.F90
    ../libflu/RFLU_DecideWrite.F90
    ../libflu/RFLU_DerivedInputValues.F90
    ../libflu/RFLU_DestroyGrid.F90
    ../libflu/RFLU_EnforceBoundsWrapper.F90
    ../libflu/RFLU_EquilibriumEulerian.F90
    ../libflu/RFLU_FinishSD.F90
    ../libflu/RFLU_GetCvLoc.F90
    ../libflu/RFLU_GetUserInput.F90
    ../libflu/RFLU_InitGlobal.F90
    ../libflu/RFLU_InitInputValues.F90
    ../libflu/RFLU_InvertMatrixSVD.F90
    ../libflu/RFLU_OpenPostInfo.F90
    ../libflu/RFLU_OpenRestartInfo.F90
    ../libflu/RFLU_PrintFlowInfo.F90
    ../libflu/RFLU_PrintFlowInfoWrapper.F90
    ../libflu/RFLU_PrintGridInfo.F90
    ../libflu/RFLU_PrintLocInfo.F90
    ../libflu/RFLU_PrintUserInput.F90
    ../libflu/RFLU_PrintWarnInfo.F90
    ../libflu/RFLU_RandomInit.F90
    ../libflu/RFLU_ReadRestartInfo.F90
    ../libflu/RFLU_ReadTbcInputFile.F90
    ../libflu/RFLU_ReadTbcSection.F90
    ../libflu/RFLU_ScalarFirst.F90
    ../libflu/RFLU_ScalarFirstPatch.F90
    ../libflu/RFLU_ScalarSecond.F90
    ../libflu/RFLU_ScalarSecondPatch.F90
    ../libflu/RFLU_ScalarCheckPositivity.F90
    ../libflu/RFLU_ScalarInitRhs.F90
    ../libflu/RFLU_ScalarViscousFluxes.F90
    ../libflu/RFLU_SetDependentVars.F90
    ../libflu/RFLU_SetGasVars.F90
    ../libflu/RFLU_SetDerivedUserInput.F90
    ../libflu/RFLU_SetModuleType.F90
    ../libflu/RFLU_SetRestartTimeFlag.F90
    ../libflu/RFLU_SetTransportVars.F90
    ../libflu/RFLU_SetVarInfo.F90
    ../libflu/RFLU_SetVarInfoWrapper.F90
    ../libflu/RFLU_SetVars.F90
    ../libflu/RFLU_SetVarsContWrapper.F90
    ../libflu/RFLU_SetVarsDiscWrapper.F90
    ../libflu/RFLU_SetVarsWrapper.F90
    ../libflu/RFLU_TestIsFirstRegion.F90
    ../libflu/RFLU_UserInput.F90
    ../libflu/RFLU_WriteRestartInfo.F90
    ../libflu/RFLU_ZeroVirtualCellVars.F90
    ../libflu/dcuhre.f)
if(RFLU_STATS)
  list(APPEND LIBFLU_SRCS
      ../libflu/RFLU_ReadStat.F90
      ../libflu/RFLU_WriteStat.F90)
endif()

set(MODFLOFLU_SRCS
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
    ../modfloflu/ModTurbulence.F90
    ../modfloflu/ModBorder.F90
    ../modfloflu/ModColoring.F90
    ../modfloflu/ModPlotting.F90)
if(RFLU_ROCSTAR)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModRocstar.F90)
endif()
if(RFLU_INRT)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesInteract.F90)
endif()
if(RFLU_PEUL)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesEulerian.F90)
endif()
if(RFLU_PLAG)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesLagrangian.F90)
endif()
if(RFLU_PERI)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesPeriodic.F90)
endif()
if(RFLU_RADI)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesRadiation.F90)
endif()
if(RFLU_SPEC)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesSpecies.F90)
endif()
if(RFLU_STATS)
  list(APPEND MODFLOFLU_SRCS
      ../modfloflu/ModInterfacesStatistics.F90
      ../modfloflu/ModStatsRoutines.F90
      ../modfloflu/ModStatsPost.F90)
endif()
if(RFLU_TURB)
  list(APPEND MODFLOFLU_SRCS ../modfloflu/ModInterfacesTurbulence.F90)
endif()

set(MODFLU_SRCS
    ../modflu/RFLU_ModAllocateMemory.F90
    ../modflu/RFLU_ModAUSMFlux.F90
    ../modflu/RFLU_ModBessel.F90
    ../modflu/RFLU_ModBFaceGradAccessList.F90
    ../modflu/RFLU_ModBilinearPatch.F90
    ../modflu/RFLU_ModBoundaryTests.F90
    ../modflu/RFLU_ModBoundLists.F90
    ../modflu/RFLU_ModBoundXvUtils.F90
    ../modflu/RFLU_ModBoundConvertCv.F90
    ../modflu/RFLU_ModCellFaceEdgeInfo.F90
    ../modflu/RFLU_ModCellMapping.F90
    ../modflu/RFLU_ModColoring.F90
    ../modflu/RFLU_ModCommLists.F90
    ../modflu/RFLU_ModConstraintUtils.F90
    ../modflu/RFLU_ModConvertCv.F90
    ../modflu/RFLU_ModCopyData.F90
    ../modflu/RFLU_ModDeallocateMemory.F90
    ../modflu/RFLU_ModDifferentiationCells.F90
    ../modflu/RFLU_ModDifferentiationFaces.F90
    ../modflu/RFLU_ModDifferentiationBFaces.F90
    ../modflu/RFLU_ModDimensionality.F90
    ../modflu/RFLU_ModDimensions.F90
    ../modflu/RFLU_ModEdgeList.F90
    ../modflu/RFLU_ModEntropyFixes.F90
    ../modflu/RFLU_ModExactFlow.F90
    ../modflu/RFLU_ModExplicitMultiStage.F90
    ../modflu/RFLU_ModFaceList.F90
    ../modflu/RFLU_ModFlowHardCode.F90
    ../modflu/RFLU_ModForcesMoments.F90
    ../modflu/RFLU_ModGatherData.F90
    ../modflu/RFLU_ModGeometry.F90
    ../modflu/RFLU_ModGeometryTools.F90
    ../modflu/RFLU_ModGrid.F90
    ../modflu/RFLU_ModGridSpeedUtils.F90
    ../modflu/RFLU_ModGridUtils.F90
    ../modflu/RFLU_ModHashTable.F90
    ../modflu/RFLU_ModHLLCFlux.F90
    ../modflu/RFLU_ModInCellTest.F90
    ../modflu/RFLU_ModInterfacesCommon.F90
    ../modflu/RFLU_ModInterfacesLibrary.F90
    ../modflu/RFLU_ModInterfacesSolver.F90
    ../modflu/RFLU_ModInterfacesUtilities.F90
    ../modflu/RFLU_ModInterpolation.F90
    ../modflu/RFLU_ModLimiters.F90
    ../modflu/RFLU_ModMPI.F90
    ../modflu/RFLU_ModNSCBC.F90
    ../modflu/RFLU_ModOctree.F90
    ../modflu/RFLU_ModOLES.F90
    ../modflu/RFLU_ModPatchCoeffs.F90
    ../modflu/RFLU_ModPatchUtils.F90
    ../modflu/RFLU_ModProbes.F90
    ../modflu/RFLU_ModReadBcInputFile.F90
    ../modflu/RFLU_ModReadWriteBcDataFile.F90
    ../modflu/RFLU_ModReadWriteFlow.F90
    ../modflu/RFLU_ModReadWriteGrid.F90
    ../modflu/RFLU_ModReadWriteGridSpeeds.F90
    ../modflu/RFLU_ModRelatedPatches.F90
    ../modflu/RFLU_ModRegionMapping.F90
    ../modflu/RFLU_ModRenumberings.F90
    ../modflu/RFLU_ModRenumberList.F90
    ../modflu/RFLU_ModRepair3D.F90
    ../modflu/RFLU_ModResidual.F90
    ../modflu/RFLU_ModRindStates.F90
    ../modflu/RFLU_ModRoeFlux.F90
    ../modflu/RFLU_ModThrustSpecImpulse.F90
    ../modflu/RFLU_ModSplitTree.F90
    ../modflu/RFLU_ModStencilsBFaces.F90
    ../modflu/RFLU_ModStencilsCells.F90
    ../modflu/RFLU_ModStencilsFaces.F90
    ../modflu/RFLU_ModStencilsVert.F90
    ../modflu/RFLU_ModStencilsUtils.F90
    ../modflu/RFLU_ModSymmetryPeriodic.F90
    ../modflu/RFLU_ModTETMESH.F90
    ../modflu/RFLU_ModTime.F90
    ../modflu/RFLU_ModTimeZoom.F90
    ../modflu/RFLU_ModTopologyUtils.F90
    ../modflu/RFLU_ModVertexLists.F90
    ../modflu/RFLU_ModViscousFlux.F90
    ../modflu/RFLU_ModWeights.F90
    ../modflu/RFLU_ModWENO.F90)
if(RFLU_ROCSTAR)
  list(APPEND MODFLU_SRCS
      ../modflu/RFLU_ModRocstarAdmin.F90
      ../modflu/RFLU_ModRocstarIO.F90
      ../modflu/RFLU_ModRocstarTools.F90
      ../modflu/RFLU_ModRocstarUtils.F90
      ../modflu/RFLU_ModInterfacesExternal.F90)
endif()
if(PETSC)
  list(APPEND MODFLU_SRCS
      ../modflu/RFLU_ModNewtonKrylov.F90
      ../modflu/RFLU_ModPETScNewtonKrylov.F90
      ../modflu/RFLU_ModPETScAdmin.F90)
endif()
if(RFLU_STATS)
  list(APPEND MODFLO_SRCS ../modflo/RFLO_ModStatsBoundaryConditions.F90)
endif()

set(TURB_SRCS
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
    ../rocturb/TURB_ModPostFlu.F90
    ../rocturb/TURB_coRansWallDistOVFlu.F90
    ../rocturb/TURB_coViscousFluxesFlu.F90
    ../rocturb/TURB_coWlmReadBcSectionFlu.F90
    ../rocturb/TURB_fluCv2Cons.F90
    ../rocturb/TURB_fluCv2Prim.F90
    ../rocturb/TURB_fluFaceVolume.F90
    ../rocturb/TURB_fluLesBLij.F90
    ../rocturb/TURB_fluLesBMij.F90
    ../rocturb/TURB_fluLesC2F.F90
    ../rocturb/TURB_rFLU_DeallocateMemory.F90
    ../rocturb/TURB_rFLU_ReadSolutionASCII.F90
    ../rocturb/TURB_rFLU_ReadSolutionBinary.F90
    ../rocturb/TURB_rFLU_WriteSolutionASCII.F90
    ../rocturb/TURB_rFLU_WriteSolutionBinary.F90)
if(RFLU_STATS)
  list(APPEND TURB_SRCS
      ../rocturb/TURB_StatCCollector.F90
      ../rocturb/TURB_StatFCollector.F90
      ../rocturb/TURB_StatMapping.F90)
endif()

set(RADI_SRCS
    ../rocrad/RADI_ModInterfaces.F90
    ../rocrad/RADI_ModParameters.F90
    ../rocrad/RADI_AllocateMemory.F90
    ../rocrad/RADI_BuildVersionString.F90
    ../rocrad/RADI_CalcEffTemp.F90
    ../rocrad/RADI_CheckParamInput.F90
    ../rocrad/RADI_DerivedInputValues.F90
    ../rocrad/RADI_DiffRadFlux.F90
    ../rocrad/RADI_DiffRadFluxPatch.F90
    ../rocrad/RADI_DiffRadIntens.F90
    ../rocrad/RADI_EmsInit.F90
    ../rocrad/RADI_ExtinctionCoef.F90
    ../rocrad/RADI_FlimConvectiveFluxes.F90
    ../rocrad/RADI_FlimDiffFlux.F90
    ../rocrad/RADI_FlimDiffFluxPatch.F90
    ../rocrad/RADI_FlimEmsInit.F90
    ../rocrad/RADI_FlimEmsUpdate.F90
    ../rocrad/RADI_FlimNumericalDissipation.F90
    ../rocrad/RADI_FlimRkInit.F90
    ../rocrad/RADI_FlimRkUpdate.F90
    ../rocrad/RADI_FlimSourceTerms.F90
    ../rocrad/RADI_FlimZeroDummyCells.F90
    ../rocrad/RADI_FluxLimiter.F90
    ../rocrad/RADI_InitInputValues.F90
    ../rocrad/RADI_InitSolution.F90
    ../rocrad/RADI_MixtSourceTermsFlim.F90
    ../rocrad/RADI_PeulSourceTermsFlim.F90
    ../rocrad/RADI_PlagSourceTermsFlim.F90
    ../rocrad/RADI_PrintUserInput.F90
    ../rocrad/RADI_ReadInputFile.F90
    ../rocrad/RADI_ReadRadiSection.F90
    ../rocrad/RADI_RkInit.F90
    ../rocrad/RADI_SolutionUpdate.F90
    ../rocrad/RADI_SourceTerms.F90
    ../rocrad/RADI_UserInput.F90)

set(PLAG_SRCS
    ../rocpart/PLAG_BuildVersionString.F90
    ../rocpart/PLAG_ModCheckVars.F90
    ../rocpart/PLAG_ModDataStruct.F90
    ../rocpart/PLAG_ModDimensions.F90
    ../rocpart/PLAG_ModEulerian.F90
    ../rocpart/PLAG_ModInjection.F90
    ../rocpart/PLAG_ModInterfaces.F90
    ../rocpart/PLAG_ModInterfacesCECells.F90
    ../rocpart/PLAG_ModParameters.F90
    ../rocpart/PLAG_ModReallocateMemory.F90
    ../rocpart/PLAG_ModRkInit.F90
    ../rocpart/PLAG_ModSurfStats.F90
    ../rocpart/PLAG_RFLU_ModComm.F90
    ../rocpart/PLAG_RFLU_ModFindCells.F90
    ../rocpart/PLAG_CalcBreakup.F90
    ../rocpart/PLAG_CalcDerivedVariables.F90
    ../rocpart/PLAG_CalcRhsPosition.F90
    ../rocpart/PLAG_CheckUserInput.F90
    ../rocpart/PLAG_DerivedInputValues.F90
    ../rocpart/PLAG_InitPatchData.F90
    ../rocpart/PLAG_InitInputValues.F90
    ../rocpart/PLAG_InjcMakeParticle.F90
    ../rocpart/PLAG_InjcSetInjection.F90
    ../rocpart/PLAG_InjcTileRKUpdate.F90
    ../rocpart/PLAG_InjcTileUpdate.F90
    ../rocpart/PLAG_InjcTileZeroRhs.F90
    ../rocpart/PLAG_IntrpMixtProperties.F90
    ../rocpart/PLAG_NonCvUpdate.F90
    ../rocpart/PLAG_PrintUserInput.F90
    ../rocpart/PLAG_ReadDisPartInitSection.F90
    ../rocpart/PLAG_ReadDisPartnContSection.F90
    ../rocpart/PLAG_ReadDisPartSection.F90
    ../rocpart/PLAG_ReadInputFile.F90
    ../rocpart/PLAG_ReadPdfFromFile.F90
    ../rocpart/PLAG_ReflectParticleData.F90
    ../rocpart/PLAG_RkInit.F90
    ../rocpart/PLAG_RkUpdateWrapper.F90
    ../rocpart/PLAG_UserInput.F90
    ../rocpart/PLAG_UpdateDataStruct.F90
    ../rocpart/PLAG_ZeroRhs.F90
    ../rocpart/PLAG_INRT_AllocMemTStep.F90
    ../rocpart/PLAG_INRT_DeallocMemTStep.F90
    ../rocpart/PLAG_RFLU_AllocMemSol.F90
    ../rocpart/PLAG_RFLU_AllocMemSolTile.F90
    ../rocpart/PLAG_RFLU_AllocMemTStep.F90
    ../rocpart/PLAG_RFLU_AllocMemTStepTile.F90
    ../rocpart/PLAG_RFLU_CorrectMixtProperties.F90
    ../rocpart/PLAG_RFLU_DeallocMemSol.F90
    ../rocpart/PLAG_RFLU_DeallocMemSolTile.F90
    ../rocpart/PLAG_RFLU_DeallocMemTStep.F90
    ../rocpart/PLAG_RFLU_DeallocMemTStepTile.F90
    ../rocpart/PLAG_RFLU_InitSolutionRandom.F90
    ../rocpart/PLAG_RFLU_InitSolutionScratch.F90
    ../rocpart/PLAG_RFLU_InitSolFromSerial.F90
    ../rocpart/PLAG_RFLU_InitSolSerialWrapper.F90
    ../rocpart/PLAG_RFLU_InitSolSerial_1D.F90
    ../rocpart/PLAG_RFLU_InitSolSerial_3D.F90
    ../rocpart/PLAG_RFLU_InjcTileCalcRhs.F90
    ../rocpart/PLAG_RFLU_InjectionDriver.F90
    ../rocpart/PLAG_RFLU_ReadSolutionASCII.F90
    ../rocpart/PLAG_RFLU_ReadSolutionBinary.F90
    ../rocpart/PLAG_RFLU_Update.F90
    ../rocpart/PLAG_RFLU_WriteSolutionASCII.F90
    ../rocpart/PLAG_RFLU_WriteSolutionBinary.F90)
if(RFLU_STATS)
  list(APPEND PLAG_SRCS ../rocpart/PLAG_StatMapping.F90)
endif()

set(PEUL_SRCS ../rocsmoke/PEUL_BuildVersionString.F90)

set(INRT_SRCS
    ../rocinteract/INRT_ModInterfaces.F90
    ../rocinteract/INRT_ModParameters.F90
    ../rocinteract/INRT_AllocateAuxillary.F90
    ../rocinteract/INRT_AugmentConSources.F90
    ../rocinteract/INRT_AugmentDisSources.F90
    ../rocinteract/INRT_BuildVersionString.F90
    ../rocinteract/INRT_BurnStatusUpdate.F90
    ../rocinteract/INRT_CalcBurning.F90
    ../rocinteract/INRT_CalcDrag.F90
    ../rocinteract/INRT_CalcHeatTransferNonBurn.F90
    ../rocinteract/INRT_CalcScouring.F90
    ../rocinteract/INRT_CheckUserInput.F90
    ../rocinteract/INRT_ComputeMaxEdges.F90
    ../rocinteract/INRT_DefineBoilingRegulation.F90
    ../rocinteract/INRT_DefineBurning.F90
    ../rocinteract/INRT_DefineDrag.F90
    ../rocinteract/INRT_DefineHeatTransferNonBurn.F90
    ../rocinteract/INRT_DefineScouring.F90
    ../rocinteract/INRT_DetermineTokens.F90
    ../rocinteract/INRT_FinishMomentumEdge.F90
    ../rocinteract/INRT_Initialize.F90
    ../rocinteract/INRT_PrintMaterialInput.F90
    ../rocinteract/INRT_PrintUserInput.F90
    ../rocinteract/INRT_ReadBoilingRegulation.F90
    ../rocinteract/INRT_ReadBurning.F90
    ../rocinteract/INRT_ReadDefaultSection.F90
    ../rocinteract/INRT_ReadDrag.F90
    ../rocinteract/INRT_ReadHeatTransferNonBurn.F90
    ../rocinteract/INRT_ReadInputFile.F90
    ../rocinteract/INRT_ReadMaterialInput.F90
    ../rocinteract/INRT_ReadScouring.F90
    ../rocinteract/INRT_SetActiveness.F90
    ../rocinteract/INRT_SetMaterial.F90
    ../rocinteract/INRT_SetParticleTemp.F90
    ../rocinteract/INRT_SetPermission.F90
    ../rocinteract/INRT_SourceTerms.F90
    ../rocinteract/INRT_UserInput.F90
    ../rocinteract/INRT_VaporEnergyConversion.F90)

set(ROCSTAR_SRCS
    ../rocstar/RFLU_CheckCouplingInput.F90
    ../rocstar/RFLU_GetBoundaryValues.F90
    ../rocstar/RFLU_GetDeformation.F90
    ../rocstar/RFLU_FlowSolverDummy.F90
    ../rocstar/RFLU_PutBoundaryValues.F90
    ../rocstar/RFLU_PutBoundaryValuesAlpha.F90
    ../rocstar/RFLU_UpdateInbuffGm.F90
    ../rocstar/Fluid_compute_integrals.F90
    ../rocstar/Fluid_finalize.F90
    ../rocstar/Fluid_preHdfOutput.F90
    ../rocstar/Fluid_postHdfOutput.F90
    ../rocstar/Rocflu_load_module.F90
    ../rocstar/Rocflu_unload_module.F90)
if(RFLU_TURB)
  list(APPEND ROCSTAR_SRCS ../rocstar/TURB_RFLU_ModRocstarAdmin.F90)
endif()
if(RFLU_STATS)
  list(APPEND ROCSTAR_SRCS
      ../rocstar/STAT_RocstarStatNaming.F90
      ../rocstar/STAT_RFLU_ModRocstarAdmin.F90)
endif()

set(INIT_SRCS
    ../utilities/rocflu/init/main.C
    ../utilities/rocflu/init/BuildVersionString.F90
    ../utilities/rocflu/init/RFLU_AllocMemSolWrapper.F90
    ../utilities/rocflu/init/RFLU_DeallocMemSolWrapper.F90
    ../utilities/rocflu/init/RFLU_InitBcDataHardCode.F90
    ../utilities/rocflu/init/RFLU_InitFlowHardCode.F90
    ../utilities/rocflu/init/RFLU_InitFlowHardCodeLim.F90
    ../utilities/rocflu/init/RFLU_InitFlowHardCodeLimWrapper.F90
    ../utilities/rocflu/init/RFLU_InitFlowHardCodeWrapper.F90
    ../utilities/rocflu/init/RFLU_InitFlowScratch.F90
    ../utilities/rocflu/init/RFLU_InitFlowScratchWrapper.F90
    ../utilities/rocflu/init/RFLU_InitFlowSerialWrapper.F90
    ../utilities/rocflu/init/RFLU_PrintHeader.F90
    ../utilities/rocflu/init/RFLU_WriteVersionString.F90
    ../utilities/rocflu/init/rfluinit.F90
    ../libflu/dcuhre.f)

set(PART_SRCS
    ../utilities/rocflu/part/main.C
    ../utilities/rocflu/part/RFLU_ModCENTAUR.F90
    ../utilities/rocflu/part/RFLU_ModCOBALT.F90
    ../utilities/rocflu/part/RFLU_ModGAMBIT.F90
    ../utilities/rocflu/part/RFLU_ModMESH3D.F90
    ../utilities/rocflu/part/RFLU_ModPartitionRegion.F90
    ../utilities/rocflu/part/RFLU_ModPartitionRegionUtils.F90
    ../utilities/rocflu/part/RFLU_ModVGRIDns.F90
    ../utilities/rocflu/part/RFLU_PrintHeader.F90
    ../utilities/rocflu/part/RFLU_ReadConvGridWrapper.F90
    ../utilities/rocflu/part/RFLU_USER_EnforcePatchCoords.F90
    ../utilities/rocflu/part/RFLU_WriteVersionString.F90
    ../utilities/rocflu/part/BuildVersionString.F90
    ../utilities/rocflu/part/rflupart.F90
    ../libflu/dcuhre.f)

set(MAP_SRCS
    ../utilities/rocflu/map/main.C
    ../utilities/rocflu/map/BuildVersionString.F90
    ../utilities/rocflu/map/RFLU_PrintHeader.F90
    ../utilities/rocflu/map/RFLU_WriteVersionString.F90
    ../utilities/rocflu/map/rflumap.F90
    ../libflu/dcuhre.f)

set(LIB_SRCS ${MODFLU_SRCS} ${MODFLOFLU_SRCS} ${LIBFLU_SRCS} ${LIBFLOFLU_SRCS} ${CORE_SRCS})

if(RFLU_ROCSTAR)
  add_definitions(-DGENX)
  list(APPEND LIB_SRCS ${ROCSTAR_SRCS})
endif()
if(RFLU_RADI)
  add_definitions(-DRADI)
  list(APPEND LIB_SRCS ${RADI_SRCS})
endif()
if(RFLU_TURB)
  add_definitions(-DTURB)
  list(APPEND LIB_SRCS ${TURB_SRCS})
endif()
if(RFLU_PLAG)
  add_definitions(-DPLAG)
  list(APPEND LIB_SRCS ${PLAG_SRCS})
endif()
if(RFLU_PEUL)
  add_definitions(-DPEUL)
  list(APPEND LIB_SRCS ${PEUL_SRCS})
endif()
if(RFLU_INRT)
  add_definitions(-DINRT)
  list(APPEND LIB_SRCS ${INRT_SRCS})
endif()
if(RFLU_STATS)
  add_definitions(-DSTATS)
endif()
if(RFLU_PERI)
  add_definitions(-DPERI)
endif()
if(RFLU_SPEC)
  add_definitions(-DSPEC)
endif()

add_library(Rocflu ${LIB_SRCS})
target_link_libraries(Rocflu Rocman Rocmop IMPACT::SITCOMF ${LAPACK_LIBRARIES})
target_include_directories(Rocflu 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rocstar>)
set_target_properties(Rocflu PROPERTIES VERSION ${ROCSTAR_VERSION}
        SOVERSION ${ROCSTAR_MAJOR_VERSION})

add_executable(rfluinit ${INIT_SRCS})
target_link_libraries(rfluinit Rocflu IMPACT::SimIN IMPACT::SimOUT)

add_executable(rflupart ${PART_SRCS})
target_link_libraries(rflupart Rocflu IMPACT::SimIN IMPACT::SimOUT ${METIS_LIB})

add_executable(rflumap ${MAP_SRCS})
target_link_libraries(rflumap Rocflu)

# Unit test programs
#read.tst: $(RTST) read.cpp
#        $(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c read.cpp
#        $(LD) $(OPTS) read.o $(RTST) -o $@ $(LFLAGS)

#adj.tst: $(ATST) adj.cpp
#        $(CXX)  $(OPTS) $(CCFLAGS) -DSTANDALONE -c adj.cpp
#        $(LD) $(OPTS) adj.o $(ATST) -o $@ $(LFLAGS)

#add_test(RunAllTests ${EXECUTABLE_OUTPUT_PATH}/testx testresults.txt)
#add_test(TestObject:ValidHandle ${TEST_RESULTS} TestObject:ValidHandle testresults.txt)
#add_test(TestObject:Exists ${TEST_RESULTS} TestObject:Exists testresults.txt)

install(TARGETS Rocflu rflupart rflumap rfluinit
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
