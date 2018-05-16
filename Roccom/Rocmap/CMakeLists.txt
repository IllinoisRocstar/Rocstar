cmake_minimum_required (VERSION 2.8)
project(Rocmap)
#ENABLE_LANGUAGE( CXX Fortran )
#include(FortranCInterface)
#FortranCInterface_HEADER(${PROJECT_SOURCE_DIR}/include/FC.h MACRO_NAMESPACE "FC_")
#FortranCInterface_VERIFY(CXX QUIET)
#set (BUILD_STATIC FALSE CACHE BOOL "Build static COM library")
#set (ENABLE_MPI TRUE CACHE BOOL "Build with MPI Support")
#set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
#set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
#mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
load_cache(../)
IF(NOT ENABLE_MPI)
  add_definitions( -DDUMMY_MPI )
ELSE()
  FIND_PACKAGE(MPI REQUIRED)
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

#INCLUDE(CTest)

#find_library(DL_LIB dl)

#set (MAPLIB_SRCS src/Rocmap.C src/Pane_boundary.C src/Pane_connectivity.C src/Pane_communicator.C 
#                 src/Dual_connectivity.C src/Simple_manifold_2.C src/Pane_ghost_connectivity.C src/KD_tree_3.C)
set (MAPLIB_SRCS src/Rocmap.C src/Pane_boundary.C src/Pane_connectivity.C src/Pane_communicator.C 
                 src/Dual_connectivity.C src/Simple_manifold_2.C src/Pane_ghost_connectivity.C)
set (UTIL_SRCS util/addpconn.C)
set (ALL_MAP_SRCS "${MAPLIB_SRCS} ${UTIL_SRCS}")

#set (TEST_SRCS test/maptest.C)
set_source_files_properties(${ALL_MAP_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC ${MPI_CXX_COMPILE_FLAGS}" )
#set(TEST_SRCS src/COMTest.C)

# rpath settings
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(include )

IF(NOT BUILD_STATIC)
  add_library(Rocmap SHARED ${MAPLIB_SRCS})
ELSE()
  add_library(Rocmap STATIC ${LIB_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
target_link_libraries(Rocmap Rocin Rocout Roccom ${MPI_CXX_LIBRARIES})

add_executable(pconntest test/pconntest.C)
target_link_libraries(pconntest Rocmap ${HDF5_LIB})
add_executable(pcommtest test/pcommtest.C)
target_link_libraries(pcommtest Rocmap ${HDF5_LIB})
add_executable(pcommpartest test/pcommpartest.C)
target_link_libraries(pcommpartest Rocmap ${HDF5_LIB})
add_executable(bordertest_hex test/bordertest_hex.C)
target_link_libraries(bordertest_hex Rocmap ${HDF5_LIB})
add_executable(bordertestg_hex test/bordertestg_hex.C)
target_link_libraries(bordertestg_hex Rocmap ${HDF5_LIB})
add_executable(bordertest_struc test/bordertest_struc.C)
target_link_libraries(bordertest_struc Rocmap ${HDF5_LIB})
add_executable(addpconn util/addpconn.C)
target_link_libraries(addpconn Rocmap ${HDF5_LIB})

IF(MPI_LINK_FLAGS)
  SET_TARGET_PROPERTIES(Rocmap pconntest pcommtest pcommpartest bordertest_hex bordertestg_hex
    bordertest_struc addpconn PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
ENDIF()

#ADD_TEST(RunAllTests ${EXECUTABLE_OUTPUT_PATH}/testx testresults.txt)
#ADD_TEST(TestObject:ValidHandle ${TEST_RESULTS} TestObject:ValidHandle testresults.txt)
#ADD_TEST(TestObject:Exists ${TEST_RESULTS} TestObject:Exists testresults.txt)

INSTALL(TARGETS Rocmap RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
