cmake_minimum_required (VERSION 2.8)
project(GridConversion)
ENABLE_LANGUAGE( CXX Fortran )
IF(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
ELSE()
  FIND_PACKAGE(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
#  IF(MPI_CXX_COMPILER)
#    set (CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#  ENDIF(MPI_CXX_COMPILER)
ENDIF()
#include(FortranCInterface)
#FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/include/FC.h MACRO_NAMESPACE "FC_")
#FortranCInterface_VERIFY(CXX QUIET)
set (BUILD_STATIC FALSE CACHE BOOL "Build static library")
set (ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

FIND_LIBRARY(IRAD_LIB IRAD)
FIND_FILE(IRAD_HDR COMM.H)
FIND_PROGRAM(RUNTEST runtest)
FIND_PROGRAM(TESTRESULT testresults)
GET_FILENAME_COMPONENT(IRAD_LIBPATH ${IRAD_LIB} PATH)
GET_FILENAME_COMPONENT(IRAD_INCPATH ${IRAD_HDR} PATH)
GET_FILENAME_COMPONENT(IRAD_EXEPATH ${RUNTEST} PATH)

file(GLOB CPP_SRCS "src/*.C")
file(GLOB F90_SRCS "src/*.f90")
#set (TESTMOD_SRCS src/TestModule.C)
#set (FTESTMOD_SRCS src/ModTestObject.f90 src/TestModule.f90)
set (ALL_SRCS "${F90_SRCS} ${CPP_SRCS}" )
set_source_files_properties(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC" )

# rpath settings
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories( include ${CMAKE_BINARY_DIR}/include ../include ${IRAD_INCPATH})


# This line runs the unit tests from the "com_test" executable
ADD_TEST(RunGridConversionUnitTests ${EXECUTABLE_OUTPUT_PATH}/gridconversion_test -o gridconversion_testresults.txt)
ADD_TEST(RunTestScripts ${RUNTEST} -l ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/tests.list -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt)
ADD_TEST(RunPlatformSpecificTests ${RUNTEST} -p ${PROJECT_SOURCE_DIR}/share/Platforms/platforms -s ${PROJECT_SOURCE_DIR} -b ${PROJECT_BINARY_DIR}/bin -o gridconversion_testresults.txt)
ADD_TEST(RunSpecificTestScript    ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Platforms/test_stub.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -o gridconversion_testresults.txt)
#Uncomment test below if you have edited and used the regression test script
#ADD_TEST(RegressionTest:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/regression.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})
ADD_TEST(TestReadingCorrectInput:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/TestReadingCorrectInput.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})
ADD_TEST(TestIncorrectInput1:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/TestIncorrectInput1.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})
ADD_TEST(TestIncorrectInput2:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/TestIncorrectInput2.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})
ADD_TEST(TestIncorrectInput3:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/TestIncorrectInput3.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})
ADD_TEST(TestIncorrectInput4:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/TestIncorrectInput4.csh -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o gridconversion_testresults.txt -a ${IRAD_EXEPATH})

# Simple test of return code of a program
ADD_TEST(ExampleProgram:Runs ${EXECUTABLE_OUTPUT_PATH}/sep ${CMAKE_BINARY_DIR}/CMakeCache.txt)

# "TESTRESULT" just checks to see if the tests worked
ADD_TEST(TestStubWorks ${TESTRESULT} TestStubWorks gridconversion_testresults.txt)
ADD_TEST(ExampleProgram:Works ${TESTRESULT} ExampleProgram:Works gridconversion_testresults.txt)
ADD_TEST(ExampleFunction:Works ${TESTRESULT} ExampleFunction:Works gridconversion_testresults.txt)
ADD_TEST(TrapezoidQuadrature:Runs ${TESTRESULT} TrapezoidQuadrature:Runs gridconversion_testresults.txt)
ADD_TEST(TrapezoidQuadrature:Accurate ${TESTRESULT} TrapezoidQuadrature:Accurate gridconversion_testresults.txt)
ADD_TEST(TrapezoidQuadrature:Order ${TESTRESULT} TrapezoidQuadrature:Order2 gridconversion_testresults.txt)
ADD_TEST(MidPointQuadrature:Runs ${TESTRESULT} MidPointQuadrature:Runs gridconversion_testresults.txt)
ADD_TEST(MidPointQuadrature:Accurate ${TESTRESULT} MidPointQuadrature:Accurate gridconversion_testresults.txt)
ADD_TEST(MidPointQuadrature:Order ${TESTRESULT} MidPointQuadrature:Order2 gridconversion_testresults.txt)
#Uncomment test below and add your Regression test name in the appropriate
#spot if you have edited and used the regression test script
#ADD_TEST(RegressionTest:Works ${TESTRESULT} YouRegressionTestNameHere gridconversion_testresults.txt)
ADD_TEST(TestReadingCorrectInput:Works ${TESTRESULT} TestReadingCorrectInput gridconversion_testresults.txt)
ADD_TEST(TestIncorrectInput1:Works ${TESTRESULT} TestIncorrectInput1 gridconversion_testresults.txt)
ADD_TEST(TestIncorrectInput2:Works ${TESTRESULT} TestIncorrectInput2 gridconversion_testresults.txt)
ADD_TEST(TestIncorrectInput3:Works ${TESTRESULT} TestIncorrectInput3 gridconversion_testresults.txt)
ADD_TEST(TestIncorrectInput4:Works ${TESTRESULT} TestIncorrectInput4 gridconversion_testresults.txt)

IF(MPI_CXX_COMPILER)
  # The parallel tests are actually run by this command
 ADD_TEST(GridConversion::RunParallelPlatformTests ${RUNTEST} -p ${PROJECT_SOURCE_DIR}/share/Platforms/parallel_platforms -s ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 3 -o gridconversion_testresults.txt)

  # The remainder use "testresults" utility to check the results of the parallel tests run above
  ADD_TEST(ParallelExample:Runs ${TESTRESULT} PEPI:Runs gridconversion_testresults.txt)
  ADD_TEST(ParallelExample:Works ${TESTRESULT} PEPI:Works gridconversion_testresults.txt)
  ADD_TEST(ParallelTrapezoidQuadrature:Runs ${TESTRESULT} ParallelTrapezoidQuadrature:Runs gridconversion_testresults.txt)
  ADD_TEST(ParallelTrapezoidQuadrature:Accurate ${TESTRESULT} ParallelTrapezoidQuadrature:Accurate gridconversion_testresults.txt)
  ADD_TEST(ParallelTrapezoidQuadrature:Order ${TESTRESULT} ParallelTrapezoidQuadrature:Order2 gridconversion_testresults.txt)
#  ADD_TEST(ParallelTrapezoidQuadrature:WeakScaling ${TESTRESULT} ParallelTrapezoidQuadrature:WeakScaling gridconversion_testresults.txt)
#  ADD_TEST(ParallelTrapezoidQuadrature:StrongScaling ${TESTRESULT} ParallelTrapezoidQuadrature:StrongScaling gridconversion_testresults.txt)
  ADD_TEST(ParallelMidPointQuadrature:Runs ${TESTRESULT} ParallelMidPointQuadrature:Runs gridconversion_testresults.txt)
  ADD_TEST(ParallelMidPointQuadrature:Accurate ${TESTRESULT} ParallelMidPointQuadrature:Accurate gridconversion_testresults.txt)
  ADD_TEST(ParallelMidPointQuadrature:Order ${TESTRESULT} ParallelMidPointQuadrature:Order2 gridconversion_testresults.txt)
#  ADD_TEST(ParallelMidPointQuadrature:WeakScaling ${TESTRESULT} ParallelMidPointQuadrature:WeakScaling gridconversion_testresults.txt)
#  ADD_TEST(ParallelMidPointQuadrature:StrongScaling ${TESTRESULT} ParallelMidPointQuadrature:StrongScaling gridconversion_testresults.txt)
ENDIF()
