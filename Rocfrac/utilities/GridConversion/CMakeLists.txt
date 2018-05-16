cmake_minimum_required(VERSION 2.8)
project(GridConversion)

# Compiled objects output paths
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
# Set up RPATH
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# 
# ========= USE IN EXISTING PROJECT ============
#
# In order to use this IllinoisRocstar Project Template
# in an existing project COPY this file from this point
# on and paste it into your project's CMakeLists.txt 
# file *before* your existing INSTALL directives.
#
INCLUDE(CTest)
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/IRAD/")
  add_subdirectory(IRAD)
  SET(IRAD_INCPATH ${CMAKE_SOURCE_DIR}/IRAD/include)
  SET(RUNTEST ${EXECUTABLE_OUTPUT_PATH}/bin/runtest)
  SET(TESTRESULT ${EXECUTABLE_OUTPUT_PATH}/bin/testresults)
  SET(IRAD_LIB IRAD)
ELSE()
  FIND_LIBRARY(IRAD_LIB IRAD)
  FIND_FILE(IRAD_HDR COMM.H)
  FIND_PROGRAM(RUNTEST runtest)
  FIND_PROGRAM(TESTRESULT testresults)
  GET_FILENAME_COMPONENT(IRAD_LIBPATH ${IRAD_LIB} PATH)
  GET_FILENAME_COMPONENT(IRAD_INCPATH ${IRAD_HDR} PATH)
  GET_FILENAME_COMPONENT(IRAD_EXEPATH ${RUNTEST} PATH)
#  link_directories(${IRAD_LIBPATH})
ENDIF()
# ADD_SUBDIRECTORY(IRAD)

IF(EXISTS ${CMAKE_SOURCE_DIR}/IMPACT)
  add_subdirectory(IMPACT)
  MESSAGE(STATUS "IMPACT will be built.")
  SET(IMPACT_INCPATH ${CMAKE_SOURCE_DIR}/IMPACT/COM/include)
  SET(IMPACT_LIB SITCOM)
  include_directories(${CMAKE_BINARY_DIR}/include)
ELSE()
  FIND_LIBRARY(IMPACT_LIB SITCOM)
  FIND_LIBRARY(SOLVERUTILS_LIB SolverUtils)
  FIND_FILE(IMPACT_HDR com.h)
  GET_FILENAME_COMPONENT(IMPACT_LIBPATH ${IMPACT_LIB} PATH)
  GET_FILENAME_COMPONENT(IMPACT_INCPATH ${IMPACT_HDR} PATH)
ENDIF()

find_package(MPI)
IF(MPI_CXX_COMPILER)
  set (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
ENDIF()

# Project-specific constructs
include_directories(include ${IRAD_INCPATH} ${IMPACT_INCPATH})

#ADD_LIBRARY(GridConversion SHARED src/ExampleSourceFile.C )
ADD_LIBRARY(GridConversion SHARED src/ExampleSourceFile.C src/GridConversion.C src/DriverProgram.C src/ReadPatranInput.C src/ReadStanfordInput.C src/HigherOrder.C src/ConnectivityMaps.C src/WriteOutput.C)
TARGET_LINK_LIBRARIES(GridConversion ${IRAD_LIB} ${IMPACT_LIB} ${SOLVERUTILS_LIB})
#INSTALL(TARGETS GridConversion LIBRARY DESTINATION lib)
ADD_EXECUTABLE(gridconversion_test src/TestGridConversion.C)
TARGET_LINK_LIBRARIES(gridconversion_test GridConversion)
ADD_EXECUTABLE(sep src/ExampleSerialProgram.C src/sep.C)
TARGET_LINK_LIBRARIES(sep GridConversion)
ADD_EXECUTABLE(make_project src/MakeProject.C)
TARGET_LINK_LIBRARIES(make_project GridConversion)
ADD_EXECUTABLE(gridConversion src/GridConversion.C src/DriverProgram.C src/ReadPatranInput.C src/ReadStanfordInput.C src/HigherOrder.C src/ConnectivityMaps.C src/WriteOutput.C)
TARGET_LINK_LIBRARIES(gridConversion GridConversion)
INSTALL(TARGETS sep make_project RUNTIME DESTINATION bin)
#INSTALL(TARGETS sep make_project gridConversion RUNTIME DESTINATION bin)

IF(MPI_CXX_COMPILER)
  ADD_DEFINITIONS(-D_GRIDCONVERSION_PARALLEL_)
  ADD_EXECUTABLE(pepi src/ExampleParallelProgram.C src/pepi.C)
  TARGET_LINK_LIBRARIES(pepi GridConversion)
  ADD_EXECUTABLE(gridconversion_parallel_test src/ParallelTestGridConversion.C)
  TARGET_LINK_LIBRARIES(gridconversion_parallel_test GridConversion)
  INSTALL(TARGETS pepi RUNTIME DESTINATION bin)
ENDIF()

# Add target to build the documentation
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ON)
IF(BUILD_DOCUMENTATION)

  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  ENDIF()
  
#  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  
  
  SET( doxyfile_in          ${CMAKE_CURRENT_SOURCE_DIR}/doc/GridConversion.cfg.in )
  SET( doxyfile             ${CMAKE_CURRENT_BINARY_DIR}/GridConversion.cfg        )
  SET( documentation        ${CMAKE_CURRENT_BINARY_DIR}/docs                 )

  CONFIGURE_FILE( ${doxyfile_in} ${doxyfile} @ONLY )

  ADD_CUSTOM_COMMAND( OUTPUT ${documentation}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
                      DEPENDS GridConversion ${doxy_extra_files}
                      COMMENT "Generating documentation")

  ADD_CUSTOM_TARGET( docs DEPENDS ${documentation} )

  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION share/GridConversion )
ENDIF()

# ========= USE IN EXISTING PROJECT ==============

ADD_SUBDIRECTORY(testing)
