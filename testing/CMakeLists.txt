cmake_minimum_required (VERSION 2.8)
project(Rocstar)
enable_language( CXX Fortran )
if(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
else()
  find_package(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
endif()
set (BUILD_STATIC FALSE CACHE BOOL "Build static library")
set (ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

find_program(RUNTEST runtest)
find_program(TESTRESULT testresults)
get_filename_component(IRAD_EXEPATH ${RUNTEST} PATH)

if(MPI_CXX_COMPILER)
  add_test(ACMRocfluRegressionTest:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/ACM_Rocflu.csh -s
      ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o rocstar_testresults.txt -a ${IRAD_EXEPATH})
  add_test(ACMRocfloRegressionTest:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/ACM_Rocflo.csh -s
      ${PROJECT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o rocstar_testresults.txt -a ${IRAD_EXEPATH})
  add_test(SSSRegressionTest:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/SSS.csh -s ${PROJECT_SOURCE_DIR}
      -b ${CMAKE_BINARY_DIR}/bin -v 2 -o rocstar_testresults.txt -a ${IRAD_EXEPATH})
  set_tests_properties(SSSRegressionTest:Runs PROPERTIES ENVIRONMENT "PATH=${IMPACT_BINPATH}:$ENV{PATH}" )
  add_test(EPRegressionTest:Runs ${RUNTEST} -f ${PROJECT_SOURCE_DIR}/share/Testing/test_scripts/EP.csh -s ${PROJECT_SOURCE_DIR}
      -b ${CMAKE_BINARY_DIR}/bin -v 2 -o rocstar_testresults.txt -a ${IRAD_EXEPATH})
  set_tests_properties(EPRegressionTest:Runs PROPERTIES ENVIRONMENT "PATH=${IMPACT_BINPATH}:$ENV{PATH}" )
  # The remainder use "testresults" utility to check the results of the parallel tests run above
  add_test(ACMRocfluRegressionTest:Works ${TESTRESULT} ACMRocfluRegressionTest rocstar_testresults.txt)
  add_test(ACMRocfloRegressionTest:Works ${TESTRESULT} ACMRocfloRegressionTest rocstar_testresults.txt)
  add_test(SSSRegressionTest:Works ${TESTRESULT} SSSRegressionTest rocstar_testresults.txt)
  add_test(EPRegressionTest:Works ${TESTRESULT} EPRegressionTest rocstar_testresults.txt)
endif()
