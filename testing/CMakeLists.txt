cmake_minimum_required(VERSION 2.8)
enable_language(CXX Fortran)

message(STATUS "ACCESSING TEST DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

set(TEST_ENV_PATH_OPTIONS "PATH=${PROJECT_BINARY_DIR}/bin:$ENV{PATH}")
set(TEST_ENV_LD_OPTIONS
    "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}/lib:\
${CMAKE_BINARY_DIR}/lib:\
${IMPACT_RUNTIME_LIBRARY_DIRS}:\
$ENV{LD_LIBRARY_PATH}")

message(STATUS "TEST_ENV_PATH_OPTIONS = ${TEST_ENV_PATH_OPTIONS}")
message(STATUS "TEST_ENV_LD_OPTIONS = ${TEST_ENV_LD_OPTIONS}")

set(BUILD_SHARED_LIBS TRUE CACHE BOOL "Build static library")
set(ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")

#adding gtest subdir with its own cmake file
add_subdirectory(lib/googletest-master)
#include gtest library. gtest_SOURCE_DIR is available from subdir addition
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

add_executable(runACMRocfluTest 
    ${PROJECT_SOURCE_DIR}/testing/src/ACM_RocfluTest.C)
target_link_libraries(runACMRocfluTest gtest gtest_main Rocflu Rocman 
    Rocmop ${MPI_CXX_LIBRARIES})

message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")

if(MPI_CXX_COMPILER)
  add_test(NAME ACM.Rocflu
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runACMRocfluTest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}
       ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/testing/
  )
else()
  message(STATUS "An MPI_CXX Compiler was not found, testing in this state is
   not possible.")
endif(MPI_CXX_COMPILER)
